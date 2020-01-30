#  GATK genomic variant calling pipeline

## New

This pipeline can be run with [Snakemake](https://snakemake.readthedocs.io/en/stable/):

[variants_GATK.snakefile](variants_GATK.snakefile) - to call variants.

[VQSR.snakefile](VQSR.snakefile) - to perform [Variant Quality Score Recalibration](https://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr).
Note, you need a good varianrt reference database.
If you do not have it, filter your variants with the [GATK VariantFiltration](https://evodify.com/gatk-in-non-model-organism/).

To execute Snakemake files locally:

```{bash eval = FALSE}
snakemake -s Snakefile --printshellcmds -j <ncores>
```

To run it on a [Slurm cluster](https://slurm.schedmd.com/overview.html):

```{bash eval = FALSE}
snakemake -s Snakefile \
    --printshellcmds \
    -j <ncores>  \
    --cluster-config cluster.yaml \
    --cluster "sbatch -A {config.account} -t {config.time} -p {config.partition} -n {config.n}"
```

## Old

The same pipeline was earlier written in Python to generate sbatch files.
These sbatch files are used on a Slurm cluster. Jobs submitted with these files are connected by dependencies when necessary.

This pipeline is described in this [blog-post](http://evodify.com/genomic-variant-calling-pipeline/)

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
