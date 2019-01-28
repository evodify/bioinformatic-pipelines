#  Generate sbatch files

Scripts to generate sbatch files and pipelines for a Slurm cluster.

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

[reads-to-VCF.py](reads-to-VCF.py) - the genomic variant calling pipeline. It is described in this [blog-post](http://evodify.com/genomic-variant-calling-pipeline/)

**DISCLAIMER:** USE THESE SCRIPTS AT YOUR OWN RISK. I MAKE NO WARRANTIES THAT THESE SCRIPTS ARE BUG-FREE, COMPLETE, AND UP-TO-DATE. I AM NOT LIABLE FOR ANY LOSSES IN CONNECTION WITH THE USE OF THESE SCRIPTS.