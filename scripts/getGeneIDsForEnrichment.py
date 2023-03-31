import sys

ids = {}

# get the ids
fh = open(sys.argv[1], 'r')
first = 0
myId = ''
newId = ''
for line in fh:
    if line[0] == ">" and first == 0:
        first = 1 
        data = line.strip().split('.')
        myId = data[-1][2:]
        i = 0
        while myId[i]  == "0":
            i += 1     
        myId = myId[i:]
        continue
    if line[0] == ">" and first == 1:
        first = 0
        newId = line.strip()[1:]
        ids[myId] = newId  
fh.close()

for myId in ids:
    print (myId, ids[myId])




