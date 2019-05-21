import time
import os
import sys

chroms = ['{}'.format(i) for i in range(1,23)]

hg19_chroms = set(['chr{}'.format(i) for i in range(1,23)])
hs37d5_chroms = set([str(i) for i in range(1,23)])
def remove_chr_from_vcf(in_vcf, out_vcf):
    with open(in_vcf, 'r') as inf, open(out_vcf, 'w') as outf:
        for line in inf:
            if line[0] == '#':
                print(line.strip(),file=outf)
                continue
            el = line.strip().split('\t')
            assert(el[0] in hg19_chroms)
            el[0] = el[0][3:]
            assert(el[0] in hs37d5_chroms)
            print("\t".join(el),file=outf)


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

if snakemake.wildcards.individual == 'NA12878' and snakemake.wildcards.build == '1000g' and snakemake.wildcards.aligner in ['blasr','blasr_filtered']:
    w_chrom = 'chr'+w_chrom
    w_ref = snakemake.input.hg19
    out_vcf = out_vcf + '.tmp'

if snakemake.wildcards.individual == 'NA12878':
    model = 'trainedModels/fullv3-pacbio-ngmlr-hg002-hg19/learningRate1e-3.epoch999'
else:
    model = 'trainedModels/fullv3-pacbio-ngmlr-hg001-hg19/learningRate1e-3.epoch999'

t1 = time.time()
cmd = '''
clairvoyante.py callVarBam \
       --chkpnt_fn {5} \
       --ref_fn {0} \
       --bam_fn {1} \
       --ctgName {2} \
       --call_fn {3} \
       --sampleName {4} \
       --threshold 0.2 \
       --minCoverage 4 \
       --threads 4
'''.format(w_ref, snakemake.input.bam, w_chrom, out_vcf, snakemake.wildcards.individual, model)
print(cmd)
s=os.system(cmd)
assert(s==0)

t2 = time.time()

if snakemake.wildcards.individual == 'NA12878' and snakemake.wildcards.build == '1000g' and snakemake.wildcards.aligner in ['blasr','blasr_filtered']:
    # remove 'chr' from reference name in vcf
    remove_chr_from_vcf(snakemake.output.vcf+'.tmp',snakemake.output.vcf)

with open(snakemake.output.runtime,'w') as outf:
    print(seconds_to_formatted_time(t2-t1),file=outf)
