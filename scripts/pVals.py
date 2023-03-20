import sys
from coalAssoc import simAndInfer

numSims = int(sys.argv[1])
popSizeArr = [0.0001,0.0002,0.0005,0.001,0.002,0.005,00.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200]

mySim = simAndInfer.SimulateILS(numSims=numSims,theta=100)
empRF = []
for N in popSizeArr:
    empRF.append(mySim.simForRFDist(N)[0])

infN = mySim.simRFDist(popSizeArr,[7.6])

# get Dist of trait vals
mySim.getTrDist(f ='../obsData/RCCdata.txt')
mySim.setTrDist()

numSims2 = int(sys.argv[2])
mySim.numSims = numSims2

# get corrected pVal dist
mySim.correctPValDist(infN[0])

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
