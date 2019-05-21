import time
import os
import sys

def seconds_to_formatted_time(ss_time):
    hh = int(ss_time / 3600)
    ss_time -= float(hh*3600)
    mm = int(ss_time / 60)
    ss_time -= float(mm*60)
    ss = int(ss_time)
    return "{}:{:02d}:{:02d}".format(hh,mm,ss)

w_chrom = snakemake.params.w_chrom
out_vcf = snakemake.output.vcf
w_ref = snakemake.params.w_ref

t1 = time.time()
cmd = 'arrow --diploid -j 16 --referenceFilename {0} --minMapQV 30 --referenceWindow {1} --algorithm quiver -o {2} {3}'.format(w_ref, w_chrom, out_vcf, snakemake.input.bam)
print(cmd)
s=os.system(cmd)
assert(s==0)

t2 = time.time()

with open(snakemake.output.runtime,'w') as outf:
    print(seconds_to_formatted_time(t2-t1),file=outf)
