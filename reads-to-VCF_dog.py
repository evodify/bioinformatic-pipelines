#! /usr/bin/env python

'''
This script generates sbatch files to submit jobs for the pipeline including read mapping, mark duplicates, BQSR,
gVCF calling and joint VCF calling.

#Example input:

/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA000787/CFA000787_S3_L001_R1_001.fastq.gz
/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA000787/CFA000787_S3_L002_R1_001.fastq.gz
/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA000787/CFA000787_S3_L003_R1_001.fastq.gz
/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA000787/CFA000787_S3_L004_R1_001.fastq.gz
/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA001614/CFA001614_S4_L001_R1_001.fastq.gz
/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA001614/CFA001614_S4_L002_R1_001.fastq.gz
/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA001614/CFA001614_S4_L003_R1_001.fastq.gz
/private/Arctiv_vs_agrarian_rawdata/170221_ST-E00215_0206_AHF35JALXX/Sample_CFA001614/CFA001614_S4_L004_R1_001.fastq.gz


#The output will consist of three sbatch files for each step of the pipeline. See the steps below.

Note! Some parts are hardcoded, e.g. the variant reference files. Modify sbatch.py file to change those parts.

#command:

$ python reads-to-VCF.py -i input.txt -r canFam3.fa -p snic2017-7-392 -n 20 -t 10-00:00:00

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import sbatch # my custom module

############################# options #############################

parser = sbatch.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the file with a list of R1.fastq file names', type=str, required=True)
parser.add_argument('-r', '--reference', help = 'name of the reference file', type=str, required=True)
parser.add_argument('-p', '--projectID', help = 'name of the project', type=str, required=True)
parser.add_argument('-n', '--ncores', help = 'maximum number of cores to request', type=str, required=True)
parser.add_argument('-t', '--time', help = 'maximum available time to request', type=str, required=True)
parser.add_argument('-a', '--aligner', help = 'aligner: BWA or stampy', type=str, required=True)

args = parser.parse_args()
if args.aligner not in ['BWA', 'stampy']:
    raise IOError('Incorrect value "%s"  in the --aligner option. Accepted values: BWA, stampy' % (args.aligner))
############################# program #############################

pathfile = open(args.input, 'r')

lanes = []
prevFullPathR1List = []
prevFullPathR1splitList = []
sampleNames = []
prevsampleName = 'NA'

for fullPathR1 in pathfile:
    fullPathR1split = fullPathR1.split("/")
    fileNameR1 = fullPathR1split[-1]
    fileNameR1L = fileNameR1.split("_")
    sampleName = fileNameR1L[0]
    laneName = fileNameR1L[2]

    if prevsampleName == 'NA' or sampleName == prevsampleName:
        lanes.append(laneName)
        prevFullPathR1List.append(fileNameR1L)
        prevFullPathR1splitList.append(fullPathR1split)
    else:
        sampleNames.append(prevsampleName)
        # sbatch script
        outputSbacth = open("reads-to-VCF_wolf.sbatch", 'w')
        outputSbacth.write("#!/bin/sh\n")
        # mapping
        jobIDs = ""
        for lane, pPath, pPathSplit in zip(lanes, prevFullPathR1List, prevFullPathR1splitList):
            outputMap = open(prevsampleName + "_" + lane + "_map.sh", 'w')
            outputMap.write("#!/bin/sh\n")
            if args.aligner == "stampy":
                sbatch.writeMapStampyJob(pPath, pPathSplit, outputMap, args.reference, prevsampleName, lane, args.ncores)
            else:
                sbatch.writeMapBWAJob(pPath, pPathSplit, outputMap, args.reference, prevsampleName, lane, args.ncores)
            outputMap.close()
            map_job = prevsampleName + "_" + lane
            sbatch.writeSbatchScript(outputSbacth, map_job, args.projectID, "20", "1-00:00:00",
                                     map_job + "_map", map_job + "_map.sh")
            jobIDs += ":$"+map_job
        lanes = [laneName]
        prevFullPathR1List = [fileNameR1L]
        prevFullPathR1splitList = [fullPathR1split]

        # merge lanes, mark duplicates, BQSR
        JobMergeMarkDuplBQSR = prevsampleName + "_mergeMarkDuplBQSR"
        sbatch.writeSbatchScript(outputSbacth, JobMergeMarkDuplBQSR, args.projectID, "1", args.time,
                                 JobMergeMarkDuplBQSR, JobMergeMarkDuplBQSR + ".sh", jobIDs)
        jobIDs = ":$"+JobMergeMarkDuplBQSR
        outputMerge = open(prevsampleName + "_mergeMarkDuplBQSR.sh", 'w')
        outputMerge.write("#!/bin/sh\n")
        sbatch.writeMergeJob(outputMerge, prevsampleName)
        sbatch.writeMarkDuplJob(outputMerge, prevsampleName)
        sbatch.writeBQSRwolfJob(outputMerge, prevsampleName, args.reference)
        sbatch.writeBQSRanalyzeCovariatesWolfJob(outputMerge, prevsampleName, args.reference)
        outputMerge.close()

        # Genotype gVCF
        JobGVCF = prevsampleName + "_gVCF"
        sbatch.writeSbatchScript(outputSbacth, JobGVCF, args.projectID, "1", args.time,
                                 JobGVCF, JobGVCF + ".sh", jobIDs)
        outputGVCF = open(prevsampleName + "_gVCF.sh", 'w')
        outputGVCF.write("#!/bin/sh\n")
        sbatch.writeHaplotypeCallerJob(outputGVCF, prevsampleName, args.ncores, args.reference)
        outputGVCF.write("\nmkdir %s\nmv %s_* %s\n" % (prevsampleName, prevsampleName, prevsampleName))
        outputGVCF.close()

        # Qualimap
        JobQualimap = prevsampleName + "_qualimap"
        sbatch.writeSbatchScript(outputSbacth, JobQualimap, args.projectID, "1", args.time,
                                 JobQualimap, JobQualimap + ".sh", jobIDs)
        outputQualimap = open(prevsampleName + "_qualimap.sh", 'w')
        outputQualimap.write("#!/bin/sh\n")
        sbatch.writeQualimaJob(outputQualimap, prevsampleName, args.ncores)
        outputQualimap.close()

    prevsampleName = sampleName
    prevFileNameR1L = fileNameR1L
    prevFullPathR1split = fullPathR1split

sampleNames.append(prevsampleName)
prevFullPathR1List.append(prevFileNameR1L)
prevFullPathR1splitList.append(prevFullPathR1split)

# mapping
jobIDs = ""
for lane, pPath, pPathSplit in zip(lanes, prevFullPathR1List, prevFullPathR1splitList):
    outputMap = open(prevsampleName + "_" + lane + "_map.sh", 'w')
    outputMap.write("#!/bin/sh\n")
    if args.aligner == "stampy":
        sbatch.writeMapStampyJob(pPath, pPathSplit, outputMap, args.reference, prevsampleName, lane, args.ncores)
    else:
        sbatch.writeMapBWAJob(pPath, pPathSplit, outputMap, args.reference, prevsampleName, lane, args.ncores)
    outputMap.close()
    map_job = prevsampleName + "_" + lane
    sbatch.writeSbatchScript(outputSbacth, map_job, args.projectID, "20", "1-00:00:00",
                             map_job + "_map", map_job + "_map.sh")
    jobIDs += ":$"+map_job
lanes = [laneName]
prevFullPathR1List = [fileNameR1L]
prevFullPathR1splitList = [fullPathR1split]

# merge lanes, mark duplicates, BQSR
JobMergeMarkDuplBQSR = prevsampleName + "_mergeMarkDuplBQSR"
sbatch.writeSbatchScript(outputSbacth, JobMergeMarkDuplBQSR, args.projectID, "1", args.time,
                         JobMergeMarkDuplBQSR, JobMergeMarkDuplBQSR + ".sh", jobIDs)
jobIDs = ":$"+JobMergeMarkDuplBQSR
outputMerge = open(prevsampleName + "_mergeMarkDuplBQSR.sh", 'w')
outputMerge.write("#!/bin/sh\n")
sbatch.writeMergeJob(outputMerge, prevsampleName)
sbatch.writeMarkDuplJob(outputMerge, prevsampleName)
sbatch.writeBQSRwolfJob(outputMerge, prevsampleName, args.reference)
sbatch.writeBQSRanalyzeCovariatesWolfJob(outputMerge, prevsampleName, args.reference)
outputMerge.close()

# Genotype gVCF
JobGVCF = prevsampleName + "_gVCF"
sbatch.writeSbatchScript(outputSbacth, JobGVCF, args.projectID, "1", args.time,
                         JobGVCF, JobGVCF + ".sh", jobIDs)
outputGVCF = open(prevsampleName + "_gVCF.sh", 'w')
outputGVCF.write("#!/bin/sh\n")
sbatch.writeHaplotypeCallerJob(outputGVCF, prevsampleName, args.ncores, args.reference)
outputGVCF.write("\nmkdir %s\nmv %s_* %s\n" % (prevsampleName, prevsampleName, prevsampleName))
outputGVCF.close()
JobGVCFid = ":$" + JobGVCF

# Qualimap
JobQualimap = prevsampleName + "_qualimap"
sbatch.writeSbatchScript(outputSbacth, JobQualimap, args.projectID, "1", args.time,
                         JobQualimap, JobQualimap + ".sh", jobIDs)
outputQualimap = open(prevsampleName + "_qualimap.sh", 'w')
outputQualimap.write("#!/bin/sh\n")
sbatch.writeQualimaJob(outputQualimap, prevsampleName, args.ncores)
outputQualimap.close()

# Join Genotyping of all gVCFs
outputGVCFall = sbatch.writeSbatchHeader(args.projectID, 1, args.time, 'GVCF', "all")
sbatch.writeGenotypeGVCFsJob(outputGVCFall, sampleNames, args.reference)
outputGVCFall.close()

pathfile.close()
