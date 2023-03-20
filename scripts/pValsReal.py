import sys
from coalAssoc import simAndInfer

mySim = simAndInfer.SimulateILS(numSims=1,theta=100)

#get real trees
mySim.getSingleCopy()

# get/set Dist of trait vals
mySim.getTrDist(f ='../obsData/RCCdata.txt')
mySim.setTrDist()

# get dist of sum stats
mySim.calcRealPVals()

# print out sum stats
i = 0
for row in mySim.real_corr:
    for gene in row:
        for thing in gene:
            print(thing, end = ' ')
        print(mySim.copyNum[i])
    i += 1

