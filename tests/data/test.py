import sys

filename = sys.argv[1]
sep = '","'
# transfos = {}
RR = {}

with open(filename, 'r') as f:
    lines = f.readlines()
    header = lines[0].split(sep)
    for line in lines[1:]:
        row = line.split(sep)
        for rule in row[10][1:-1].split(','):
            RR[rule] = row[2]
        # transfos[row[1]] = row[2]

# filename = sys.argv[2]
# rules = {}
# with open(filename, 'r') as f:
#     lines = f.readlines()
#     header = lines[0].split(sep)
#     for line in lines[1:]:
#         row = line.split(sep)
#         rules[row[2]] = transfos[row[1][:-2]]


for rule, smiles in RR.items():
    print(rule, smiles)