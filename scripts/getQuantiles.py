import sys

fh = open(sys.argv[1],'r')
pV = float(sys.argv[2])
lim = int(sys.argv[3])

# get cc dist for each species
data = {}

i  = 0
j = 0
for line in fh:
    l = line.strip().split()
    if l[4] in data:
        data[l[4]].append([float(l[2]),j])
    else:
        data[l[4]] = [[float(l[2]),j]]
    i += 1
    if i == 8:
        i = 0
        j += 1 

# sort each dist
sortData = {}

for spec in data:
    sortOrd = sorted(data[spec], key=lambda x: x[0])
    sortData[spec] = sortOrd

# get p-vals for each i
pVals = {}
for spec in sortData:
    i = 0
    for loc in sortData[spec]:
        if loc[1] in pVals:
            pVals[loc[1]].append(i/100000.)
        else:
            pVals[loc[1]] = [i/100000.]

        i += 1

tot = 0
for i in pVals:
    j = 0
    for val in pVals[i]:
        if val < pV:
            j += 1
    if j >= lim:
        tot += 1    

print(tot)
