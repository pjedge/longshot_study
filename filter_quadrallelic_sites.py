import sys

for line in sys.stdin:
    if line[0] == '#':
        print(line,end='')
    else:
        el = line.strip().split('\t')
        alts = el[4]
        if len(alts.split(',')) <= 2:
            print(line,end='')
            
