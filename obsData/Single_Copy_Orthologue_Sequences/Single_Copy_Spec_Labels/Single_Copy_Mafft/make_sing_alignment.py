
# 
import sys
import os

seqs = {}

dir_list = os.listdir(".")
for f in dir_list:
    if not (f[-2:] == "fa" and f[:2] == "OG"):
        continue
    fh = open(f, 'r')
    spec = ''
    for line in fh:
        if line[0] == ">":
            spec = line[1:].strip()
            if spec not in seqs:
                seqs[spec] = ''
            continue
        seqs[spec] += line.strip()
    fh.close()

for spec in seqs:
    print('>'+spec)
    i = 0
    while i < len(seqs[spec]):
        print(seqs[spec][i:(i+80)])
        i += 80


