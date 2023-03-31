

import sys

# get ids
fh = open(sys.argv[1], 'r')

ids = {}

for line in fh:
    data = line.strip().split()
    ids[data[0]] = data[1]
fh.close()


# get foreground
fh = open(sys.argv[2], 'r')

idsF = []
for line in fh:
    data = line.strip()
    idsF.append(data)

for ID in idsF:
    if ID in ids:
        print(ids[ID])
