import sys
from coalAssoc import simAndInfer

mySim = simAndInfer.SimulateILS()

for edge in mySim.sp_tree.preorder_edge_iter():
    if edge.length is not None:
        edge.length = edge.length/0.0699

print(mySim.sp_tree.as_string(schema="newick",))

#mySim.setPopSize()
#mySim.initRandomTree()
#mySim.makeTaxonMap()
#mySim.simGeneTree()
#mySim.simCausalGene()
#mySim.writeSpecTree()
#mySim.makeOutputHeader()
#mySim.runNullSims()
