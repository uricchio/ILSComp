import sys
from coalAssoc import simAndInfer

numSims = int(sys.argv[1])
popSizeArr = [0.0001,0.0002,0.0005,0.001,0.002,0.005,00.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200]

mySim = simAndInfer.SimulateILS(numSims=numSims,theta=100)
empRF = []
for N in popSizeArr:
    empRF.append(mySim.simForRFDist(N)[0])

infN = mySim.simRFDist(popSizeArr,empRF)

for val in range(len(infN)):
    print (infN[val], popSizeArr[val])
