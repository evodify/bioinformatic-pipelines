#! /usr/bin/env python

'''
This script generates sbatch files to submit jobs for dogs read mapping.

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

args = parser.parse_args()

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
        # mapping
        outputMap = sbatch.writeSbatchHeader(args.projectID, args.ncores, "1-00:00:00", prevsampleName, "map")
        for lane, pPath, pPathSplit in zip(lanes, prevFullPathR1List, prevFullPathR1splitList):
            sbatch.writeMapBWAJob(pPath, pPathSplit, outputMap, args.reference, prevsampleName, lane, args.ncores)
        sbatch.writeNextSbath(outputMap, prevsampleName, "mergeMarkDuplBQSR")
        outputMap.close()
        lanes = [laneName]
        prevFullPathR1List = [fileNameR1L]
        prevFullPathR1splitList = [fullPathR1split]

        # merge lanes, mark duplicates, BQSR
        outputMerge = sbatch.writeSbatchHeader(args.projectID, 1, args.time, prevsampleName, "mergeMarkDuplBQSR")
        sbatch.writeMergeJob(outputMerge, prevsampleName)
        sbatch.writeMarkDuplJob(outputMerge, prevsampleName)
        sbatch.writeBQSRJob(outputMerge, prevsampleName, args.reference)
        sbatch.writeNextSbath(outputMerge, prevsampleName, "gVCF")
        #sbatch.writeNextSbath(outputMerge, prevsampleName, "qualimap")
        sbatch.writeBQSRanalyzeCovariatesJob(outputMerge, prevsampleName, args.reference)
        outputMerge.close()

        # Qualimap
        #outputQualimap = sbatch.writeSbatchHeader(args.projectID, args.ncores, args.time, prevsampleName, "qualimap")
        #sbatch.writeQualimaJob(outputQualimap, prevsampleName, args.ncores)
        #outputQualimap.close()

        # Genotype gVCF
        outputGVCF = sbatch.writeSbatchHeader(args.projectID, args.ncores, args.time, prevsampleName, "gVCF")
        sbatch.writeHaplotypeCallerJob(outputGVCF, prevsampleName, args.ncores, args.reference)
        outputGVCF.write("\nmkdir %s\nmv %s_* %s\n" % (prevsampleName, prevsampleName, prevsampleName))
        outputGVCF.close()

    prevsampleName = sampleName
    prevFileNameR1L = fileNameR1L
    prevFullPathR1split = fullPathR1split

sampleNames.append(prevsampleName)
prevFullPathR1List.append(prevFileNameR1L)
prevFullPathR1splitList.append(prevFullPathR1split)

# mapping
outputMap = sbatch.writeSbatchHeader(args.projectID, args.ncores, "1-00:00:00", prevsampleName, "map")
for lane, pPath, pPathSplit in zip(lanes, prevFullPathR1List, prevFullPathR1splitList):
    sbatch.writeMapBWAJob(pPath, pPathSplit, outputMap, args.reference, prevsampleName, lane, args.ncores)
sbatch.writeNextSbath(outputMap, prevsampleName, "mergeMarkDuplBQSR")
outputMap.close()

# merge lanes, mark duplicates, BQSR
outputMerge = sbatch.writeSbatchHeader(args.projectID, 1, args.time, prevsampleName, "mergeMarkDuplBQSR")
sbatch.writeMergeJob(outputMerge, prevsampleName)
sbatch.writeMarkDuplJob(outputMerge, prevsampleName)
sbatch.writeBQSRJob(outputMerge, prevsampleName, args.reference)
sbatch.writeNextSbath(outputMerge, prevsampleName, "gVCF")
#sbatch.writeNextSbath(outputMerge, prevsampleName, "qualimap")
sbatch.writeBQSRanalyzeCovariatesJob(outputMerge, prevsampleName, args.reference)
outputMerge.close()

# Qualimap
#outputQualimap = sbatch.writeSbatchHeader(args.projectID, args.ncores, args.time, prevsampleName, "qualimap")
#sbatch.writeQualimaJob(outputQualimap, prevsampleName, args.ncores)
#outputQualimap.close()

# Genotype gVCF
outputGVCF = sbatch.writeSbatchHeader(args.projectID, args.ncores, args.time, prevsampleName, "gVCF")
sbatch.writeHaplotypeCallerJob(outputGVCF, prevsampleName, args.ncores, args.reference)
outputGVCF.write("\nmkdir %s\nmv %s_* %s\n" % (prevsampleName, prevsampleName, prevsampleName))
outputGVCF.close()

# Join Genotyping of all gVCFs
outputGVCFall = sbatch.writeSbatchHeader(args.projectID, 1, "1-00:00:00", 'GVCF', "all")
sbatch.writeGenotypeGVCFsJob(outputGVCFall, sampleNames, args.reference)
outputGVCFall.close()

pathfile.close()
