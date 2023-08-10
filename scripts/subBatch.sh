#!/bin/bash	
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 4 # cores requested
#SBATCH --mem=8000 # memory in Mb
#SBATCH --array=1-1 # memory in Mb
#SBATCH -t 10:00:00 # time
#SBATCH -o outfile # send stdout to outfile
#SBATCH -e errfile # send stderr to errfile

nSpec=$1
nGene=$2

for i in {1..200}; do python simGroupGenes.py 200 100000 $nSpec $nGene >> ../simData/simGroup/highDiscP.nLin${nSpec}.nG${nGene}.txt; done
for i in {1..200}; do python simGroupGenes.py 200 0.001 $nSpec $nGene >> ../simData/simGroup/lowDiscP.nLin${nSpec}.nG${nGene}.txt; done
for i in {1..200}; do python simGroupGenes.py 200 1 $nSpec $nGene >> ../simData/simGroup/midDiscP.nLin${nSpec}.nG${nGene}.txt; done
