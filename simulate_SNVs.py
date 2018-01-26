import copy
import random
import os
import pysam
import itertools
import numpy as np
from numpy.random import choice

VALID_CHROMS = set(['chr{}'.format(c) for c in range(1,23)]+['chrX','chrY','chrM'])

# estimate prior probability of genotypes using strategy described here:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/
# "prior probability of each genotype"

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION!
# Selects a variant genotype, assuming the site already is a SNV.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def create_phased_genotype_selector():

    alleles = ['A','C','G','T']
    genotypes = list(itertools.combinations_with_replacement(alleles,2))
    het_snp_rate = 0.0005
    hom_snp_rate = 0.001
    diploid_genotype_priors = dict()
    haploid_genotype_priors = dict()
    transition = {'A':'G','G':'A','T':'C','C':'T'}

    for allele in alleles:

        # priors on haploid alleles
        haploid_genotype_priors[allele] = dict()
        haploid_genotype_priors[allele][allele] = 1 - hom_snp_rate
        haploid_genotype_priors[allele][transition[allele]] = hom_snp_rate / 6 * 4
        for transversion in alleles:
            if transversion in haploid_genotype_priors[allele]:
                continue
            haploid_genotype_priors[allele][transversion] =  hom_snp_rate / 6

        diploid_genotype_priors[allele] = []
        for G in genotypes:
            g1,g2 = G
            # probability of homozygous reference is the probability of neither het or hom SNP
            if g1 == g2 and g1 == allele:
                diploid_genotype_priors[allele].append(0)
            elif g1 == g2 and g1 != allele:
                # transitions are 4 times as likely as transversions
                if g1 == transition[allele]:
                    diploid_genotype_priors[allele].append(het_snp_rate / 6 * 4)
                else:
                    diploid_genotype_priors[allele].append(het_snp_rate / 6)
            else: # else it's the product of the haploid priors
                    diploid_genotype_priors[allele].append(haploid_genotype_priors[allele][g1] * haploid_genotype_priors[allele][g2])

        # remove the option of selecting homozygous reference
        total = sum(diploid_genotype_priors[allele])
        for i in range(len(genotypes)):
            diploid_genotype_priors[allele][i] /= total
        # make sure everything sums to 1
        diploid_genotype_priors[allele][-1] = 1.0 - sum(diploid_genotype_priors[allele][:-1])


    g_ixs = list(range(len(genotypes)))
    def phased_genotype_selector(ref_allele):
        g_ix = choice(g_ixs, 1, p=diploid_genotype_priors[ref_allele])[0]
        g = genotypes[g_ix]
        if random.random() > 0.5:
            return g
        else:
            return (g[1],g[0])

    return phased_genotype_selector


def simulate_SNV_VCF(hg19_fasta, output_vcf):

    phased_genotype_selector = create_phased_genotype_selector()

    with pysam.FastaFile(hg19_fasta) as fasta, open(output_vcf,'w') as outv:

        for chrom in fasta.references:
            size = fasta.get_reference_length(chrom)

            pos = 0
            while pos < size:

                pos += np.random.geometric(0.0015, size=None)

                ref_allele = fasta.fetch(chrom,pos,pos+1).upper()

                if ref_allele in ['A','C','G','T']:
                    genotype = phased_genotype_selector(ref_allele)
                else:
                    continue

                #hap1_seq.append(genotype[0])
                #hap2_seq.append(genotype[1])
                if genotype[0] == genotype[1] and genotype[0] != ref_allele:
                    var_str = genotype[0]
                    genotype_str = '1|1'
                elif genotype[1] == ref_allele and genotype[0] != ref_allele:
                    var_str = genotype[0]
                    genotype_str = '1|0'
                elif genotype[0] == ref_allele and genotype[1] != ref_allele:
                    var_str = genotype[1]
                    genotype_str = '0|1'
                elif genotype[0] != ref_allele and genotype[1] != ref_allele and genotype[0] != genotype[1]:
                    # triallelic
                    var_str = genotype[0] + ',' + genotype[1]
                    genotype_str = '1|2'
                else:
                    print("INVALID GENOTYPE ENCOUNTERED")
                    exit(1)

                assert(chrom in VALID_CHROMS)

                if genotype != (ref_allele,ref_allele) and 'N' not in genotype:

                    el = [None]*10
                    el[0] = chrom
                    el[1] = pos + 1
                    el[2] = '.'
                    el[3] = ref_allele
                    el[4] = var_str
                    el[5] = 100
                    el[6] = 'PASS'
                    el[7] = 'None'
                    el[8] = 'GT'
                    el[9] = genotype_str
                    line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(*el)
                    print(line,file=outv)


if __name__ == '__main__':
    generate_fasta()
