import sys


# command:
# python makeQQplotdata.py ~/projects/ILSSims/ILSsims/obsData/corrPvalDist.txt ~/projects/ILSSims/ILSsims/obsData/realPvalDist.txt  > ~/projects/ILSSims/ILSsims/obsData/corrPvalDistReal.txt > ~/projects/ILSSims/ILSsims/obsData/QQplotData.txt

# read in null corr coeffs
# store array for each species and sort to get pVals
null = {}
fh = open(sys.argv[1], 'r')
for line in fh:
    data = line.strip().split()
    if data[0] in null:
        null[data[0]].append(float(data[1]))
    else:
        null[data[0]] = [float(data[1])]
fh.close() 

for spec in null:
    null[spec] = sorted(null[spec])

# read in obs arr and store
obs = {}
fh = open(sys.argv[2], 'r')

for line in fh:
    data = line.strip().split()
    if data[2] in obs:
        obs[data[2]].append([float(data[0]),int(data[3]), float(data[1])])
    else:
        obs[data[2]] = [[float(data[0]),int(data[3]), float(data[1])]]
fh.close()

for spec in obs:
    obs[spec] = sorted(obs[spec],key=lambda x: x[0])

#  Walk through the theoretical array to get pVals
for spec in obs:
    posObs = 0    
    posNull = 0
    for item in obs[spec]:
        while posNull < len(null[spec]) and item[0] > null[spec][posNull]:
            posNull += 1  
        print(spec, item[0], item[1], max([1-(1+posNull)/(len(null[spec])),1/len(null[spec])]), max([1-(1+posObs)/(len(obs[spec])+1),1/len(obs[spec])]), end =' ')
        if float(item[0] < 0):
            print (1 - item[2]/2)
        else:
            print(item[2]/2)
    
        posObs += 1


