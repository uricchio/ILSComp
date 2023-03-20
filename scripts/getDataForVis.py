import sys
from coalAssoc import simAndInfer

mySim = simAndInfer.SimulateILS(numSims=1,theta=100)

#get real trees
mySim.getSingleCopy(myMin=5719,myMax=5719)

# get/set Dist of trait vals
mySim.getTrDist(f ='../obsData/RCCdata.txt')
mySim.setTrDist()

# get dist of sum stats
mySim.getIndVals()


