import sys
from coalAssoc import simAndInfer

param = sys.argv[1]
simDistFile =  '../obsData/pValsNull/pValsNull.'+param+'.txt'
realDistFile = '../obsData/realPvalDist.txt'

mySim = simAndInfer.SimulateILS()
mySim.calcSumStats(param,simDistFile,realDistFile)

#mySim.setPopSize()
#mySim.initRandomTree()
#mySim.makeTaxonMap()
#mySim.simGeneTree()
#mySim.simCausalGene()
#mySim.writeSpecTree()
#mySim.makeOutputHeader()
#mySim.runNullSims()
