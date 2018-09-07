#!/usr/bin/python2

'''
This is a python module to create sbatch files for job submission on a cluster
'''

############################# modules #############################

import argparse, sys # for input options

############################# classes  ############################

class CommandLineParser(argparse.ArgumentParser):
    '''To check the input arguments '''
    def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

############################# functions ###########################

def writeSbatchHeader(projectID, cores, time, sampleName, sbatchName):
    '''Creates an output file and writes the header to the file'''
    outputFile = open(sampleName + "_" + sbatchName + ".sbatch", 'w')
    outputFile.write("#!/bin/bash\n#SBATCH -A %s\n#SBATCH -p core\n#SBATCH -n %s\n#SBATCH -t %s\n#SBATCH -J %s_%s\n"
                 "#SBATCH -o %s_%s.out\n#SBATCH -e %s_%s.err\n"
                 % (projectID, cores, time, sampleName, sbatchName, sampleName, sbatchName, sampleName, sbatchName))
    return outputFile

def writeSbatchScript(outputFile, jobID, projectID, cores, time, jobName, scriptName, dependency="none"):
    '''Writes a script line to submit sbatch job and takes into account dependencies if necessary'''
    if dependency=="none":
        outputFile.write("\n%s=$(sbatch -A %s -p core -n %s -t %s -J %s -e %s.err -o %s.out %s | cut -d \" \" -f 4)\n" %
                         (jobID, projectID, cores, time, jobName, jobName, jobName, scriptName))
    else:
        outputFile.write("\n%s=$(sbatch --dependency=afterok%s -A %s -p core -n %s -t %s -J %s -e %s.err -o %s.out %s "
                         "| cut -d \" \" -f 4)\n" %
                         (jobID, dependency, projectID, cores, time, jobName, jobName, jobName, scriptName))
    outputFile.write("echo \"%s under the ID $%s has been submitted\"\n" % (scriptName, jobID))

def writeMapBWAJob(fileNameR1list, filePathR1list, outputFile, reference, sample, laneID, cores, ram):
    '''Writes the BWA alignment job to the file'''
    # make R1, R2 file path first:
    filePathR1 = '/'.join(str(e) for e in filePathR1list) # R1
    #R2:
    fileNameR1list[3] = "R2"
    fileNameR2 = '_'.join(str(e) for e in fileNameR1list)
    filePathR1list[-1] = fileNameR2
    filePathR2 = '/'.join(str(e) for e in filePathR1list)
    # remove trailing \n
    filePathR1 = filePathR1.rstrip()
    filePathR2 = filePathR2.rstrip()
    # write:
    outputFile.write("\nbwa mem -t %s -M -R '@RG\\tID:%s_%s\\tPL:illumina\\tLB:%s\\tSM:%s' %s %s %s "
                     "| samtools sort -O bam -@ %s -m %sG > %s_%s.bam\n"
                     % (cores, sample, laneID, sample, sample, reference, filePathR1, filePathR2, cores, ram, sample, laneID))

def writeMapStampyJob(fileNameR1list, filePathR1list, outputFile, divergence, reference, sample, laneID, cores, ram, tmp):
    '''Writes the stampy alignment job to the file'''
    # make R1, R2 file path first:
    filePathR1 = '/'.join(str(e) for e in filePathR1list) # R1
    #R2:
    fileNameR1list[3] = "R2"
    fileNameR2 = '_'.join(str(e) for e in fileNameR1list)
    filePathR1list[-1] = fileNameR2
    filePathR2 = '/'.join(str(e) for e in filePathR1list)
    # remove trailing \n
    filePathR1 = filePathR1.rstrip()
    filePathR2 = filePathR2.rstrip()
    referenceStampy = reference.split(".")[0]
    outputFile.write("\nstampy.py -g %s -h %s --substitutionrate=%s -t %s -o %s/%s_%s_stampy.sam -M %s %s\n"
                     "samtools sort -O bam -@ %s -m %sG %s/%s_%s_stampy.sam > %s/%s_%s.bam\n"
                     "\njava -Xmx%sG -jar picard.jar AddOrReplaceReadGroups \\\n"
                     "I=%s/%s_%s.bam \\\n"
                     "O=%s_%s.bam \\\n"
                     "RGID=%s_%s \\\n"
                     "RGLB=%s \\\n"
                     "RGPL=illumina \\\n"
                     "RGPU=%s_%s \\\n"
                     "RGSM=%s\n"
                     % (referenceStampy, referenceStampy, divergence, cores, tmp, sample, laneID, filePathR1, filePathR2,
                        cores, ram,  tmp, sample, laneID, tmp, sample, laneID,
                        ram, tmp, sample, laneID, sample, laneID, sample, laneID, sample, sample, laneID, sample))

def writeNextSbath(outputFile, sample, jobName):
    '''Writes an sbatch line to submit the next job when this job is finished.'''
    outputFile.write("\nsbatch %s_%s.sbatch\n" % (sample, jobName))


def writeMergeJob(outputFile, sample, tmp):
    '''Writes the merge job lines to the file'''
    outputFile.write("\nls %s_*_L00?.bam | sed 's/  / /g;s/ /\\n/g' > %s_bam.list\n"
                     "\nsamtools merge -b %s_bam.list %s/%s_merged.bam\n"
                     "\nxargs -a %s_bam.list rm\n" %
                     (sample, sample, sample, tmp, sample, sample))


def writeMarkDuplJob(outputFile, sample, ram, tmp):
    '''Writes the MarkDuplicates job lines to the file'''
    outputFile.write("\njava -Xmx%sG -Djava.io.tmpdir=%s -jar picard.jar "
                     "MarkDuplicates \\\nVALIDATION_STRINGENCY=LENIENT \\\nMETRICS_FILE=%s/%s_merged_markDupl_metrix.txt \\\n"
                     "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=15000 \\\nINPUT=%s/%s_merged.bam \\\n"
                     "OUTPUT=%s/%s_merged_markDupl.bam\n"
                     "\nrm %s/%s_merged.bam\n" %
                     (ram, tmp, tmp, sample, tmp, sample, tmp, sample,  tmp, sample))

def writeBQSRJob(outputFile, sample, reference, ram, tmp, known_sites):
    '''Writes the BaseRecalibrator job lines to the file'''
    outputFile.write("\n# Generate the first pass BQSR table file\n"
                     "gatk --java-options \"-Xmx%sG\"  BaseRecalibrator \\\n"
                     "-R %s \\\n"
                     "-I %s/%s_merged_markDupl.bam \\\n" % (ram, reference, tmp, sample))
    known_sitesP = known_sites.split(",")
    for f in known_sitesP:
        outputFile.write("--known-sites %s \\\n" % f)
    outputFile.write("-O %s_merged_markDupl_BQSR.table\n"
                     "\n# Apply BQSR\n"
                     "gatk --java-options \"-Xmx%sG\"  ApplyBQSR \\\n"
                     "-R %s \\\n"
                     "-I %s/%s_merged_markDupl.bam \\\n"
                     "-bqsr %s_merged_markDupl_BQSR.table \\\n"
                     "-O %s_merged_markDupl_BQSR.bam\n"
                     "\nrm %s/%s_merged_markDupl.bam\n"
                     % (sample, ram, reference, tmp, sample, sample, sample, tmp, sample))

def writeBQSRanalyzeCovariatesJob(outputFile, sample, reference, ram, known_sites):
    '''Writes the BaseRecalibrator results check job lines to the file'''
    outputFile.write("\n# Generate the second pass BQSR table file\n"
                     "gatk --java-options \"-Xmx%sG\"  BaseRecalibrator \\\n"
                     "-R %s \\\n"
                     "-I %s_merged_markDupl_BQSR.bam \\\n"
                     % (ram, reference, sample))
    known_sitesP = known_sites.split(",")
    for f in known_sitesP:
        outputFile.write("--known-sites %s \\\n" % f)
    outputFile.write("-O %s_merged_markDupl_BQSR2.table\n"
                     "\n# Plot the recalibration results\n"
                     "gatk --java-options \"-Xmx%sG\"  AnalyzeCovariates \\\n"
                     "-before %s_merged_markDupl_BQSR.table \\\n"
                     "-after %s_merged_markDupl_BQSR2.table \\\n"
                     "-plots %s_merged_markDupl_BQSR.pdf\n"
                     % (sample, ram, sample, sample, sample))

def writeQualimapJob(outputFile, sample, cores, ram):
    '''Writes the samtools index job line to the file'''
    outputFile.write("\nqualimap bamqc -nt %s --java-mem-size=%sG -bam %s_merged_markDupl_BQSR.bam -outdir %s_qualimap\n"
                     % (cores, ram, sample, sample))

def writeHaplotypeCallerJob(outputFile, sample, cores, ram, reference):
    '''Writes the HaplotypeCaller in GVCF mode job lines to the file'''
    outputFile.write("\ngatk --java-options \"-Xmx%sG\" HaplotypeCaller \\\n"
                     "-R %s \\\n"
                     "-ERC GVCF \\\n"
                     "-I %s_merged_markDupl_BQSR.bam \\\n"
                     "-O %s_merged_markDupl_BQSR.g.vcf.gz\n"
                     "\ntabix %s_merged_markDupl_BQSR.g.vcf.gz\n"
                     % (ram, reference, sample, sample, sample))


def writeGenotypeGVCFsJob(outputFile, samples, ram, reference):
    '''Writes the GenotypeGVCFs job lines to the file'''
    outputFile.write("\ngatk --java-options \"-Xmx%sg -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" CombineGVCFs \\\n"
                     "-R %s \\\n" % (ram, reference))
    for sample in samples:
        outputFile.write("-V %s/%s_merged_markDupl_BQSR.g.vcf.gz \\\n" % (sample, sample))
    outputFile.write("-O GVCF_merged_markDupl_BQSR.g.vcf.gz\n"
                     "\ntabix GVCF_merged_markDupl_BQSR.g.vcf.gz\n")

    outputFile.write("\ngatk --java-options \"-Xmx%sg -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" GenotypeGVCFs \\\n"
                     "-R %s \\\n"
                     "-V GVCF_merged_markDupl_BQSR.g.vcf.gz \\\n"
                     "-O GVCF_merged_markDupl_BQSR.vcf.gz\n"
                     "\ntabix GVCF_merged_markDupl_BQSR.vcf.gz\n" % (ram, reference))
