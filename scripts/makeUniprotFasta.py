import sys

fh = open(sys.argv[1], 'r')

for line in fh:
    data = line.strip().split('\t')
    if len(data) < 3:
        continue
    gene_names = data[2]
    names = gene_names.split()
    final_name = '>' 
    for name in names:
        if name[0:2] == "PCH":
            final_name += name[4:]
            break
        elif name[0] == "P":
            final_name += name
            break
    if final_name == ">":
        continue
    print(final_name)
    print(data[1])
   




