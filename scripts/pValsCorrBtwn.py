import sys
from coalAssoc import simAndInfer

popSizeFile  = '../obsData/popSize.txt'
popSize = 0
fh = open(popSizeFile, 'r')
for line in fh:
    popSize = float(line.strip())
    break

mySim = simAndInfer.SimulateILS(theta=100)

# get Dist of trait vals
mySim.getTrDist(f ='../obsData/RCCdata.txt')
mySim.setTrDist()

numSims = int(sys.argv[1])
mySim.numSims = numSims

# get corrected pVal dist
mySim.correctPValDist(popSize,save=False)

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
