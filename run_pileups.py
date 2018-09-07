import sys
import os

# first argument: bam file
# second argument: positions file

bam_file = sys.argv[1]
with open(sys.argv[2],'r') as inf:
    for line in inf:
        el = line.strip().split('\t')
        assert(len(el) == 2)
        chrom,pos = el
        os.system("samtools mpileup -s {0} --region {1}:{2}-{2} 2>/dev/null".format(bam_file,chrom,pos))
