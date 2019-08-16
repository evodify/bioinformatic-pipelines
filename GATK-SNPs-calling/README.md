#  GATK genomic variant calling pipeline (sbatch)

The pipeline is described in this [blog-post](http://evodify.com/genomic-variant-calling-pipeline/)

These scripts generate sbatch files for a Slurm cluster. Jobs submitted with these sbatch files will be connected by dependencies when necessary.

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

[reads-to-VCF.py](reads-to-VCF.py) - scrip that generates sbatch files for the genomic variant calling pipeline.

**DISCLAIMER:** USE THESE SCRIPTS AT YOUR OWN RISK. I MAKE NO WARRANTIES THAT THESE SCRIPTS ARE BUG-FREE, COMPLETE, AND UP-TO-DATE. I AM NOT LIABLE FOR ANY LOSSES IN CONNECTION WITH THE USE OF THESE SCRIPTS.
