import sys
from coalAssoc import simAndInfer

mySim = simAndInfer.SimulateILS(numSims=1,theta=100)

#get real trees
perc,pos =  mySim.getCorrPVals(nullDistFile='../obsData/corrPvalDistRuby.txt',trueCorrCoeffFile='../obsData/realPvalDistRuby.txt')

for spec in perc:
    for i in range(len(perc[spec])):
        print(spec, pos[spec][i], perc[spec][i])
