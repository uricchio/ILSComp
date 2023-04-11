import sys
from coalAssoc import simAndInfer

mySim = simAndInfer.SimulateILS()
mySim. getSingleCopy()
m_sd = mySim.calcRFDist()

# mean and sd of RF dist of GT compared to spec trees
print (m_sd[0], m_sd[1]) 

