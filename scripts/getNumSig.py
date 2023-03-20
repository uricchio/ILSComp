

import sys

fh = open(sys.argv[1], 'r')

num = 0
numSig = 0
for line in fh:
    data = line.strip().split()
    if float(data[3]) < 0.05:
        numSig +=1 
    num += 1
    if num == 8:
        print(numSig)
        num = 0
        numSig = 0
