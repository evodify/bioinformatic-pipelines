#!/bin/bash -l
#SBATCH -A snic2019-x-xxx
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J sbatchSnakefile
#SBATCH -e sbatchSnakefile.err
#SBATCH -o sbatchSnakefile.out

snakemake -s Snakefile -j 100 --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} -t {cluster.time} -p {cluster.partition} -n {cluster.n}"
