import sys

filename = sys.argv[1]
sep = '\t'
# transfos = {}
RR = {}

with open(filename, 'r') as f:
    lines = f.readlines()
    header = lines[0].split(sep)
    for line in lines[1:]:
        row = line.split(sep)
        print(row[0], row[7]+'>>'+row[9])
        # for rule in row[10][1:-1].split(','):
        #     RR[rule] = row[2]


for rule, smiles in RR.items():
    print(rule, smiles)