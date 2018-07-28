import sys

chroms = ['chr{}'.format(x) for x in range(1,23)]
chroms += ['{}'.format(x) for x in range(1,23)]

assert(sys.argv[-1] in ['filter_bed_chroms.py','--remove_chr'])

for line in sys.stdin:
    el = line.strip().split()
    if el[0] not in chroms:
        continue
    if sys.argv[-1] == '--remove_chr' and line[:3] == 'chr':
        print(line[3:],end='')
    else:
        print(line,end='')
