#!/bin/bash     
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 4 # cores requested
#SBATCH --mem=8000 # memory in Mb
#SBATCH --array=1-1 # memory in Mb
#SBATCH -t 10:00:00 # time
#SBATCH -o outfile # send stdout to outfile
#SBATCH -e errfile # send stderr to errfile

./run_disc_sims.sh
