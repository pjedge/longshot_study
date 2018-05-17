import sys

chroms = ['chr{}'.format(x) for x in range(1,23)] + ['chrX','chrY']

for line in sys.stdin:
    el = line.strip().split()
    if el[0] not in chroms:
        continue
    print(line[3:],end='')
