#! /usr/bin/env python
'''
This script generates sbatch files to submit jobs for the pipeline including
read mapping, mark duplicates, BQSR, gVCF calling and joint VCF calling.
BQSR step can be skipped if -k is not specified.


#Example input:

Sample_CFA000787/CFA000787_S3_L001_R1_001.fastq.gz
Sample_CFA000787/CFA000787_S3_L002_R1_001.fastq.gz
Sample_CFA000787/CFA000787_S3_L003_R1_001.fastq.gz
Sample_CFA000787/CFA000787_S3_L004_R1_001.fastq.gz
Sample_CFA001614/CFA001614_S4_L001_R1_001.fastq.gz
Sample_CFA001614/CFA001614_S4_L002_R1_001.fastq.gz
Sample_CFA001614/CFA001614_S4_L003_R1_001.fastq.gz
Sample_CFA001614/CFA001614_S4_L004_R1_001.fastq.gz


where CFA00078 and CFA001614 are sample names;
L001 - L004 are different lanes;
R1 read 1 of paired-end reads. Do not specify names with R2.

The output will consist of three sbatch files for each step of the pipeline.
See the steps below.

Note! This script requires the module sbatch.py.
So, download and place it in the same directory.

You can find more instructions on how to use this script at:
http://evodify.com/genomic-variant-calling-pipeline/

#command:

python reads-to-VCF.py \
-i R1.txt \
-r canFam3.fa \
-p snic2017 \
-n 20 \
-t 5-00:00:00 \
-l \$SNIC_TMP \
-a BWA \
-m 4 \
-k dogs.557.publicSamples.ann.chrAll.PASS.vcf.gz,00-All_chrAll.vcf.gz

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

import sbatch  # my custom module

parser = sbatch.CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help = 'name of the file with a list of R1.fastq file names',
    type = str,
    required = True)
parser.add_argument(
    '-r',
    '--reference',
    help = 'name of the reference file',
    type = str,
    required = True)
parser.add_argument(
    '-p',
    '--projectID',
    help = 'name of the project',
    type = str,
    required = True)
parser.add_argument(
    '-n',
    '--ncores',
    help = 'maximum number of cores to request (default 2)',
    type = int,
    required = False,
    default = 2)
parser.add_argument(
    '-t',
    '--time',
    help = 'maximum available time to request (default 1-00:00:00)',
    type = str,
    required = False,
    default = "1-00:00:00")
parser.add_argument(
    '-a',
    '--aligner',
    help = 'aligner: BWA or stampy (default BWA)',
    type = str,
    required = False,
    default = "BWA")
parser.add_argument(
    '-d',
    '--divergence',
    help = 'mean divergence from the references for Stampy (default 0.001)',
    type = float,
    required = False,
    default = 0.001)
parser.add_argument(
    '-m',
    '--memory',
    help = 'Maximum available RAM per core in GB (default 4)',
    type = int,
    required = False,
    default = 4)
parser.add_argument(
    '-k',
    '--knowSites',
    help  =  'comma separated list of known variants for BQSR',
    type = str,
    required = False,
    default = None)
parser.add_argument(
    '-l',
    '--tmpDirectory',
    help = 'path to the temporary directory',
    type = str,
    required = True)

args = parser.parse_args()

if args.aligner not in ['BWA', 'stampy']:
    raise IOError(
        'Incorrect value "%s" in the --aligner option.'
        'Accepted values: BWA, stampy' % (args.aligner))

def loopProcessSbatch(lanes, prevFullPathR1List, prevFullPathR1splitList,
                      aligner, divergence, reference,
                      ncores, memory, tmpDirectory,
                      time, projectID, knowSites):
    jobIDs = ""
    for lane, pPath, pPathSplit in zip(lanes, prevFullPathR1List,
                                        prevFullPathR1splitList):
        outputMap = open(prevsampleName + "_" + lane + "_map.sh", 'w')
        outputMap.write("#!/bin/sh\n")
        if args.aligner == "stampy":
            sbatch.writeMapStampyJob(pPath, pPathSplit, outputMap,
                                     divergence, reference,
                                     prevsampleName, lane, ncores,
                                     memory, tmpDirectory)
        else:
            sbatch.writeMapBWAJob(pPath, pPathSplit, outputMap,
                                  reference, prevsampleName, lane,
                                  ncores, memory)
        outputMap.close()
        map_job = prevsampleName + "_" + lane
        sbatch.writeSbatchScript(outputSbatch, map_job, projectID,
                                 ncores, time, map_job +
                                 "_map", map_job + "_map.sh")
        jobIDs += ":$" + map_job
    lanes = [laneName]
    prevFullPathR1List = [fileNameR1L]
    prevFullPathR1splitList = [fullPathR1split]

    # merge lanes, mark duplicates, BQSR
    nameTail = sbatch.ifBQSR(knowSites)
    JobMergeMarkDuplBQSR = prevsampleName + "_" + nameTail

    sbatch.writeSbatchScript(
        outputSbatch, JobMergeMarkDuplBQSR, projectID,
        "1", time, JobMergeMarkDuplBQSR,
        JobMergeMarkDuplBQSR + ".sh", jobIDs)

    jobIDs = ":$" + JobMergeMarkDuplBQSR
    
    outputMerge = open(prevsampleName + "_" + nameTail + ".sh", 'w')
    outputMerge.write("#!/bin/sh\n")
    sbatch.writeMergeJob(outputMerge, prevsampleName, tmpDirectory)
    sbatch.writeMarkDuplJob(outputMerge, prevsampleName, memory,
                            tmpDirectory)
    if args.knowSites!=None:
        sbatch.writeBQSRJob(outputMerge, prevsampleName, args.reference,
                            memory, tmpDirectory, knowSites)
        sbatch.writeBQSRanalyzeCovariatesJob(outputMerge, prevsampleName,
                                             reference, memory, knowSites)
    outputMerge.close()

    # Genotype gVCF
    JobGVCF = prevsampleName + "_gVCF"
    sbatch.writeSbatchScript(outputSbatch, JobGVCF, projectID,
                                "1", time, JobGVCF,
                                JobGVCF + ".sh", jobIDs)
    outputGVCF = open(prevsampleName + "_gVCF.sh", 'w')
    outputGVCF.write("#!/bin/sh\n")
    sbatch.writeHaplotypeCallerJob(outputGVCF, prevsampleName, ncores,
                                   memory, reference, knowSites)
    outputGVCF.write("\nmkdir %s\n"
                        "mv %s_* %s\n" %
                        (prevsampleName,
                        prevsampleName, prevsampleName))
    outputGVCF.close()

    # Qualimap
    JobQualimap = prevsampleName + "_qualimap"
    sbatch.writeSbatchScript(outputSbatch, JobQualimap, projectID,
                                ncores, time, JobQualimap,
                                JobQualimap + ".sh", jobIDs)
    outputQualimap = open(prevsampleName + "_qualimap.sh", 'w')
    outputQualimap.write("#!/bin/sh\n")
    sbatch.writeQualimapJob(outputQualimap, prevsampleName, ncores,
                            memory, knowSites)
    outputQualimap.close()
    return [[laneName], [fileNameR1L], [fullPathR1split]]

# program

pathfile = open(args.input, 'r')

lanes = []
prevFullPathR1List = []
prevFullPathR1splitList = []
sampleNames = []
prevsampleName = 'NA'

# sbatch script
outputSbatch = open("reads-to-VCF.sbatch", 'w')
outputSbatch.write("#!/bin/sh\n")

laneName = 'L1'

# loop through the input fatsq list
for fullPathR1 in pathfile:
    fullPathR1split = fullPathR1.split("/")
    fileNameR1 = fullPathR1split[-1]
    fileNameR1L = fileNameR1.split("_")
    sampleName = fileNameR1L[0]
    laneName = '_'.join(str(e) for e in fileNameR1L[1:3])
   

    # check if the sample name is the same but the lane number if different
    if prevsampleName == 'NA' or sampleName == prevsampleName:
        lanes.append(laneName)
        prevFullPathR1List.append(fileNameR1L)
        prevFullPathR1splitList.append(fullPathR1split)
    else:
        sampleNames.append(prevsampleName)
        outLPS = loopProcessSbatch(
                    lanes, prevFullPathR1List, prevFullPathR1splitList,
                    args.aligner, args.divergence, args.reference,
                    args.ncores, args.memory, args.tmpDirectory,
                    args.time, args.projectID, args.knowSites)
        lanes, prevFullPathR1List, prevFullPathR1splitList = outLPS

    prevsampleName = sampleName
    prevFileNameR1L = fileNameR1L
    prevFullPathR1split = fullPathR1split

outLPS = loopProcessSbatch(lanes, prevFullPathR1List, prevFullPathR1splitList,
                           args.aligner, args.divergence, args.reference,
                           args.ncores, args.memory, args.tmpDirectory,
                           args.time, args.projectID, args.knowSites)

# Join Genotyping of all gVCFs
sampleNames.append(prevsampleName)
outputGVCFall = sbatch.writeSbatchHeader(args.projectID, 1, args.time,
                                         'GVCF', "all")
sbatch.writeGenotypeGVCFsJob(outputGVCFall, sampleNames, args.memory,
                             args.reference, args.knowSites)
outputGVCFall.close()

pathfile.close()
