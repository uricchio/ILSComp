import sys
from coalAssoc import simAndInfer

numSims = int(sys.argv[1])
popSize = float(sys.argv[2])
theta = 10/popSize # float(sys.argv[3])

mySim = simAndInfer.SimulateILS(popSize = popSize, theta=theta,numSims=numSims)
mySim.simTraitAndNull()

#mySim.setPopSize()
#mySim.initRandomTree()
#mySim.makeTaxonMap()
#mySim.simGeneTree()
#mySim.simCausalGene()
#mySim.writeSpecTree()
#mySim.makeOutputHeader()
#mySim.runNullSims()
