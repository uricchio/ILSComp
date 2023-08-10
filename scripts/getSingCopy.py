import sys
from coalAssoc import simAndInfer

mySim = simAndInfer.SimulateILS(numSims=1,theta=100)

#get real trees
mySim.getSingleCopy()
for tr in mySim.myTrees:
    mySim.outputCladogram(tr)

