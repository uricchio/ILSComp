import sys
from coalAssoc import simAndInfer

popSize = float(sys.argv[1])
mySim = simAndInfer.SimulateILS(numSims=10000,theta=100)

# get Dist of trait vals
mySim.getTrDist(f ='../obsData/RCCdata.txt')
mySim.setTrDist()

numSims2 = int(sys.argv[2])
mySim.numSims = numSims2

# get corrected pVal dist
mySim.correctPValDist(popSize)

for sp in mySim.nullCorr:
    for val in mySim.nullCorr[sp]:
        print(sp[1:-1], val)

#mySim.setPopSize()
#mySim.initRandomTree()
#mySim.makeTaxonMap()
#mySim.simGeneTree()
#mySim.simCausalGene()
#mySim.writeSpecTree()
#mySim.makeOutputHeader()
#mySim.runNullSims()
