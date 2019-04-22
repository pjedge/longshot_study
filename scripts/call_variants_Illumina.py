import time
import os
import sys

ref_file = {'1000g':'data/genomes/hs37d5.fa', 'hg38':'data/genomes/hg38.fa'}
chroms = ['{}'.format(i) for i in range(1,23)]

def seconds_to_formatted_time(ss_time):
    hh = int(ss_time / 3600)
    ss_time -= float(hh*3600)
    mm = int(ss_time / 60)
    ss_time -= float(mm*60)
    ss = int(ss_time)
    return "{}:{:02d}:{:02d}".format(hh,mm,ss)

# this function takes in a chromosome name in 1..22,X and a "genome build" name
# and appends a 'chr' prefix to the chrom name if the genome build is hg38
def chr_prefix(chrom, build):
    assert(build in ['1000g','hg38'])
    assert(chrom in chroms)
    if build == 'hg38':
        return 'chr'+chrom
    else:
        return chrom

w_chrom = chr_prefix(snakemake.wildcards.chrom, snakemake.wildcards.build)
w_ref = ref_file[snakemake.wildcards.build]
t1 = time.time()
cmd = '''
freebayes -f {} \
--standard-filters \
--region {} \
 --genotype-qualities \
 {} \
  > {}
'''.format(w_ref, w_chrom, snakemake.input.bam, snakemake.output.vcf)

s=os.system(cmd)
assert(s==0)

t2 = time.time()
with open(snakemake.output.runtime,'w') as outf:
    print(seconds_to_formatted_time(t2-t1),file=outf)
