#
import sys
import os 

fList = os.listdir(sys.argv[1])

keep = 0

for f in fList:
    if 'fa' not in f:
        continue
    fh = open(os.path.join(sys.argv[1],f), 'r')
    for line in fh:
    
        if ">" in line and keep == 1:
            keep = 0
            break

        if '53503-Pc' in line:
            print(line.strip()+'.'+f[:-3])
            keep = 1
            continue
        if keep == 1:
            print(line.strip())
    fh.close()
