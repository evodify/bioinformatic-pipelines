#  Generate sbatch files

Scripts to generate many sbatch for job submitting on Slurm cluster and pipelinining.

The content of `script.sbatch` looks like this:
```
#!/bin/bash
#SBATCH -A project-name
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J job-name
#SBATCH -o output.out
#SBATCH -e errors.err

command-to-execute
```

[sbatch.py](sbatch.py) is a custom module that contains all necessary functions. It is loaded as `import sbatch `.

[reads-to-VCF_dog.py](reads-to-VCF_dog.py) - the genomic variant calling pipeline:
*mapping lanes -> merge BAM files from all lanes -> mark duplicates -> perfrom BQSR -> check mapping quality -> genotype in GVCF mode-> joint genotyping of all gVCF files*.
It mostly follows [the GATK best practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145). See the details in the file. This particular pipeline was used to genotype dogs.
