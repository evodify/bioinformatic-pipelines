#!/usr/bin/python2

'''
This is a python module to create sbacth files for job submission on a cluster
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
    outputFile = open(sbatchName + "_" + sampleName + ".sbatch", 'w')
    outputFile.write("#!/bin/bash\n#SBATCH -A %s\n#SBATCH -p core\n#SBATCH -n %s\n#SBATCH -t %s\n#SBATCH -J %s_%s\n"
                 "#SBATCH -o %s_%s.out\n#SBATCH -e %s_%s.err\n"
                 % (projectID, cores, time, sbatchName, sampleName, sbatchName, sampleName, sbatchName, sampleName))
    return outputFile


def writeMapJob(fileNameR1list, filePathR1list, outputFile, reference, sample, laneID, cores):
    '''Writes the job lines to the file'''
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
                     "| samtools sort -O bam -@ %s > %s_%s.bam\n"
                     % (cores, sample, laneID, sample, sample, reference, filePathR1, filePathR2, cores, sample, laneID))


def writeNextSbath(outputFile, sample, jobName):
    '''Writes an sbatch line to submit the next job when this job is finished.'''
    outputFile.write("\nsbatch %s_%s.sbatch\n" % (jobName, sample))


def writeMergeJob(outputFile, sample):
    '''Writes the merge job lines to the file'''
    outputFile.write("\nls %s_*.bam | sed 's/  / /g;s/ /\\n/g' > %s_bam.list\n"
                     "\nsamtools merge -b %s_bam.list $SNIC_TMP/%s_merged.bam\n"
                     "\nxargs -a %s_bam.list rm\nrm %s_bam.list\n" %
                     (sample, sample, sample, sample, sample, sample))


def writeMarkDuplJob(outputFile, sample):
    '''Writes the MarkDuplicates job lines to the file'''
    outputFile.write("\njava -Xmx6g -Djava.io.tmpdir=$SNIC_TMP -jar /sw/apps/bioinfo/picard/2.10.3/rackham/picard.jar "
                     "MarkDuplicates \\\nVALIDATION_STRINGENCY=LENIENT \\\nMETRICS_FILE=%s_merged_markDupl_metrix.txt \\\n"
                     "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=15000 \\\nINPUT=$SNIC_TMP/%s_merged.bam \\\n"
                     "OUTPUT=$SNIC_TMP/%s_merged_markDupl.bam\n"
                     "\nsamtools index $SNIC_TMP/%s_merged_markDupl.bam\n" %
                     (sample, sample, sample, sample))

def writeBQSRJob(outputFile, sample, reference):
    '''Writes the BaseRecalibrator job lines to the file'''
    outputFile.write("\n# Generate the first pass BQSR table file\n"
                     "gatk --java-options \"-Xmx6G\"  BaseRecalibrator \\\n"
                     "-R %s \\\n"
                     "-I $SNIC_TMP/%s_merged_markDupl.bam \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/Axelsson_2013_SNPs_canfam3.vcf \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/Axelsson_2013_indels_canfam3.vcf \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/160DG.99.9.recalibrated_variants.vcf.gz \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/00-All_chrAll.vcf \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/dogs.557.publicSamples.ann.chrAll.PASS.vcf \\\n"
                     "-O %s_merged_markDupl_BQSR.table\n"
                     "\n# Apply BQSR\n"
                     "gatk --java-options \"-Xmx6G\"  ApplyBQSR \\\n"
                     "-R %s \\\n"
                     "-I $SNIC_TMP/%s_merged_markDupl.bam \\\n"
                     "-bqsr %s_merged_markDupl_BQSR.table \\\n"
                     "-O %s_merged_markDupl_BQSR.bam\n"
                     % (reference, sample, sample, reference, sample, sample, sample))


def writeBQSRanalyzeCovariatesJob(outputFile, sample, reference):
    '''Writes the BaseRecalibrator results check job lines to the file'''
    outputFile.write("\n# Generate the second pass BQSR table file\n"
                     "gatk --java-options \"-Xmx6G\"  BaseRecalibrator \\\n"
                     "-R %s \\\n"
                     "-I %s_merged_markDupl_BQSR.bam \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/Axelsson_2013_SNPs_canfam3.vcf \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/Axelsson_2013_indels_canfam3.vcf \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/160DG.99.9.recalibrated_variants.vcf.gz \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/00-All_chrAll.vcf \\\n"
                     "--known-sites /proj/uppstore2017236/b2013119/private/dmytro/BQSR-reference/dogs.557.publicSamples.ann.chrAll.PASS.vcf \\\n"
                     "-O %s_merged_markDupl_BQSR2.table\n"
                     "\n# Plot the recalibration results\n"
                     "gatk --java-options \"-Xmx6G\"  AnalyzeCovariates \\\n"
                     "-before %s_merged_markDupl_BQSR.table \\\n"
                     "-after %s_merged_markDupl_BQSR2.table \\\n"
                     "-plots %s_merged_markDupl_BQSR.pdf\n"
                     % (reference, sample, sample, sample, sample, sample))

def writeQualimaJob(outputFile, sample, cores):
    '''Writes the samtools index job line to the file'''
    outputFile.write("\nqualimap bamqc -nt %s --java-mem-size=6G -bam %s_merged_markDupl_BQSR.bam -outdir %s_qualimap\n"
                     % (cores, sample, sample))

def writeHaplotypeCallerJob(outputFile, sample, cores, reference):
    '''Writes the HaplotypeCaller in GVCF mode job lines to the file'''
    outputFile.write("\ngatk --java-options \"-Xmx6G\" HaplotypeCaller \\\n"
                     "-native-pair-hmm-threads %s \\\n"
                     "-R %s \\\n"
                     "-ERC GVCF \\\n"
                     "-I %s_merged_markDupl_BQSR.bam \\\n"
                     "-O %s_merged_markDupl_BQSR.g.vcf\n"
                     % (cores, reference, sample, sample))


def writeGenotypeGVCFsJob(outputFile, samples, reference):
    '''Writes the GenotypeGVCFs job lines to the file'''
    outputFile.write("\ngatk --java-options \"-Xmx6G\" GenotypeGVCFs \\\n"
                     "-R %s \\\n" % reference)
    for sample in samples:
        outputFile.write("-V %s_merged_markDupl_BQSR.g.vcf \\\n" % sample)
    outputFile.write("-O GVCF_merged_markDupl_BQSR.vcf\n")
