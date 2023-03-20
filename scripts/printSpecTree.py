import sys
from coalAssoc import simAndInfer

mySim = simAndInfer.SimulateILS()

print(mySim.sp_tree.as_string(schema="nexus",))

#mySim.setPopSize()
#mySim.initRandomTree()
#mySim.makeTaxonMap()
#mySim.simGeneTree()
#mySim.simCausalGene()
#mySim.writeSpecTree()
#mySim.makeOutputHeader()
#mySim.runNullSims()
