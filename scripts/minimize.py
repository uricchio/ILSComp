import sys
from coalAssoc import simAndInfer

popSizeStatFile =  '../obsData/pValsNull/sumDataDiffs.txt'

mySim = simAndInfer.SimulateILS()
mySim.interpPopSize(popSizeStatFile)

