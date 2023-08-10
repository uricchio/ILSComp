for nLin in {10,20,50,100,200,500}; do for nGenes in {1,2,5,10,20}; do sbatch subBatch.sh $nLin $nGenes; done; done
