import sys
from coalAssoc import simAndInfer

numSims = int(sys.argv[1])
popSize = float(sys.argv[2])
nLin = int(sys.argv[3])
nGene = int(sys.argv[4])
theta = 10/popSize # float(sys.argv[3])

mySim = simAndInfer.SimulateILS(popSize = popSize, theta=theta,numSims=numSims)
mySim.simPowerSpecTree(nLin=nLin,nGene=nGene)

#mySim.setPopSize()
#mySim.initRandomTree()
#mySim.makeTaxonMap()
#mySim.simGeneTree()
#mySim.simCausalGene()
#mySim.writeSpecTree()
#mySim.makeOutputHeader()
#mySim.runNullSims()
