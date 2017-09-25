import re
import sys
from bisect import bisect_left
import pysam
import pickle
from matplotlib import pyplot as plt


bases = {'A','T','G','C'}


def has_homopolymer(mystr, N):
    curr_letter = mystr[0]
    count = 1

    for letter in mystr[1:]:
        if letter == curr_letter:
            count += 1
        else:
            count = 1
            curr_letter = letter

        if count >= N:
            return True

    return False

# make a dictionary from VCF where the key is genome coordinates and value is genotype
def make_vcf_dict(vcf_file,prepend_chr=False, SNPs_only=True, remove_homopolymers=False, hg19=None, min_qual=0, min_gq=0, min_dp=0):

    vcf_dict = dict()
    dp_pat = re.compile("DP=(\d+)")

    if remove_homopolymers:
        hg19_fasta = pysam.FastaFile(hg19)

    with open(vcf_file,'r') as inf:
        for line in inf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]

            if prepend_chr:
                vcf_chrom = 'chr' + vcf_chrom

            vcf_pos = int(vcf_line[1])

            ref_allele = str.upper(vcf_line[3])
            variant_allele = str.upper(vcf_line[4])
            third_allele = ''

            if remove_homopolymers and has_homopolymer(str.upper(hg19_fasta.fetch(vcf_chrom, vcf_pos-5-1, vcf_pos+5-1)), 3):
                continue

            fields = vcf_line[7]

            if float(vcf_line[5]) < min_qual:
                continue

            if min_dp > 0:
                depth = int(float(re.findall(dp_pat, fields)[0]))
                if depth < min_dp:
                    continue

            if vcf_line[8][0:2] != 'GT':
                continue
            genotype = vcf_line[9][0:3]
            if genotype in ['0/0','0|0']:
                continue

            if min_gq > 0:
                gq_ix = None
                for i,label in enumerate(vcf_line[8].split(':')):
                    if label == 'GQ':
                        gq_ix = i
                        break

                gq = float(vcf_line[9].split(':')[gq_ix])

                if gq < min_gq:
                    continue

            if not genotype[1] in ['/','|'] and genotype[0] in ['0','1','2'] and genotype[2] in ['0','1','2']:
                continue

            if ',' in variant_allele:
                splt = variant_allele.split(',')
                if len(splt) != 2:
                    continue
                variant_allele, third_allele = splt

            if SNPs_only and len(ref_allele) > 1 or len(variant_allele) > 1 or len(third_allele) > 1:
                continue

            # check that all 3 alleles are valid strings of bases
            assert(set(list(ref_allele + variant_allele + third_allele)) <= bases)

            alleles = [None,None]

            if genotype[0] == '0':
                alleles[0] = ref_allele
            elif genotype[0] == '1':
                alleles[0] = variant_allele
            elif genotype[0] == '2':
                assert(third_allele != '') # check that we really had a third allele
                alleles[0] = third_allele
            else:
                print(line.strip())
                print("ERROR: Invalid genotype in VCF file")
                exit(1)

            if genotype[2] == '0':
                alleles[1] = ref_allele
            elif genotype[2] == '1':
                alleles[1] = variant_allele
            elif genotype[2] == '2':
                assert(third_allele != '') # check that we really had a third allele
                alleles[1] = third_allele
            else:
                print(line.strip())
                print("ERROR: Invalid genotype in VCF file")
                exit(1)

            g = tuple(sorted(alleles))

            vcf_dict[(vcf_chrom,vcf_pos)] = (ref_allele, g)

    return vcf_dict

# make a dictionary from HapCUT2 output where the key is genome coordinates and value is genotype
def make_hap_dict(hap_variant_calls, cutoff=-1):

    # read in the calls we made with PacBio reads
    hap_dict = dict()

    with open(hap_variant_calls,'r') as inf:
        for line in inf:
            if 'BLOCK' in line:
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue
            if el[1] == '-':
                continue
            ref_allele = el[5]
            var_allele = el[6]

            hap_alleles = {el[1],el[2]}
            chrom = el[3]
            pos   = int(el[4])

            a1, a2 = None,None
            if hap_alleles == {'0'}:
                a1, a2 = (ref_allele, ref_allele)
                continue # for now, we ignore reference calls
            elif hap_alleles == {'0','1'}:
                a1, a2 = sorted((ref_allele, var_allele))
            elif hap_alleles == {'1'}:
                a1, a2 = (var_allele, var_allele)
            else:
                print(line.strip())
                print("ERROR: Invalid genotype in haplotype file")
                exit(1)

            g = (a1, a2) # ref_allele, var_allele, a1, a2) --- a1 and a2 are unordered

            if float(el[10]) >= cutoff:
                hap_dict[(chrom,pos)] = (ref_allele, g)

    return hap_dict

# simple class to represent bed intervals
class bed_file:

    regions = dict()
    region_starts = dict()

    def __init__(self, bed_fn, prepend_chr=False, chrom_filter=None):

        def parse_bedfile(input_file):

            boundaries = []
            with open(input_file,'r') as inf:
                for line in inf:

                    if len(line) < 3:
                        continue

                    el = line.strip().split('\t')

                    chrom = el[0]

                    if prepend_chr:
                        chrom = 'chr'+chrom

                    start = int(el[1])
                    stop  = int(el[2])

                    boundaries.append((chrom, start, stop))

            return boundaries

        bedlist = parse_bedfile(bed_fn)
        if chrom_filter == None:
            chroms = set([x[0] for x in bedlist])
        else:
            chroms = [chrom_filter]

        for c in chroms:
            self.regions[c] = [(start,stop) for (chrom,start,stop) in bedlist if chrom == c]
            self.region_starts[c] = [start for (chrom,start,stop) in bedlist if chrom == c]

    def covered(self,chrom,pos):
        ix = bisect_left(self.region_starts[chrom],pos)
        start,stop = self.regions[chrom][ix - 1]
        if pos >= start and pos < stop:
            return True
        else:
            return False

    def coverage(self, chrom):
        total = 0
        for chrom, boundlist in self.regions.items():
            total += sum([stop - start for (start,stop) in boundlist])

        return total


# compute a confusion matrix for a set of test variant calls vs set of reference variant calls (GIAB).
def compare_test_ref(test, ref, bed_filter, chrom_filter,hg19_handle=None):

    TP = 0
    FP = 0
    FN = 0
    FP_errors = set()
    misaligned = 0

    for (chrom, pos), (ref_allele,(a1,a2)) in test.items():

        if not bed_filter.covered(chrom,pos):
            continue

        if (chrom,pos) in ref:
            if ref[(chrom,pos)] == (ref_allele,(a1,a2)):
                TP += 1 # correctly called variant
            else:
                #print("FP at {} {}".format(chrom,pos))
                FP_errors.add((chrom,pos))
                FP += 1 # called variant but not the right one

        else: # missed this variant
            FP_errors.add((chrom, pos))
            FP += 1

    for (chrom, pos), (ref_allele,(a1,a2)) in ref.items():

        if not bed_filter.covered(chrom,pos):
            continue

        if (chrom,pos) not in test:
            FN += 1

    TN = bed_filter.coverage(chrom_filter) - (TP + FP + FN)
    return (TP, TN, FP, FN), FP_errors

def count_accuracy_old(realigned_vcf, illumina_vcf, giab_ref_vcf, bed_filter_file, chrom_filter, hg19_fasta, report):

    # files with variants we want to test the accuracy of
    quals = [0,10,20,30,40,50,60,70,80,90,100,200,300,400]


    # read the variant calls for the reference dataset
    giab_ref = make_vcf_dict(giab_ref_vcf)

    bed_filter = bed_file(bed_filter_file,chrom_filter=chrom_filter)

    # compute confusion matrix for each set of calls
    realigned_precisions = []
    realigned_recalls    = []
    illumina_precisions  = []
    illumina_recalls     = []

    #for qual in quals:
    qual = 30

    print(qual)
    #print("TESTING PACBIO REALIGNMENT:")
    realigned_test = make_vcf_dict(realigned_vcf, min_gq=50)

    realigned_res, realigned_FP_errors = compare_test_ref(realigned_test,giab_ref,bed_filter, chrom_filter)
    #if ((realigned_res[0] + realigned_res[2]) == 0 or (realigned_res[0]+realigned_res[3]) == 0):
    #    break
    realigned_precisions.append(realigned_res[0] / (realigned_res[0] + realigned_res[2]))
    realigned_recalls.append(realigned_res[0]/(realigned_res[0]+realigned_res[3]))

    #print("TESTING STANDARD ILLUMINA:")
    illumina_test = make_vcf_dict(illumina_vcf,min_qual=qual)
    illumina_res, illumina_FP_errors = compare_test_ref(illumina_test,giab_ref,bed_filter, chrom_filter)
    #if ((illumina_res[0] + illumina_res[2]) == 0 or (illumina_res[0]+illumina_res[3]) == 0):
    #    break
    illumina_precisions.append(illumina_res[0] / (illumina_res[0] + illumina_res[2]))
    illumina_recalls.append(illumina_res[0]/(illumina_res[0]+illumina_res[3]))

    #print("")
    #if qual == quals[0]:
    with open(report,'w') as outf:
        for current_file in [outf,sys.stdout]:
            print("pacbio realignment\tTP:{}\tTN:{}\tFP:{}\tFN:{}".format(*realigned_res),file=current_file)
            print("precision: {} recall: {}".format(realigned_res[0]/(realigned_res[0]+realigned_res[2]), realigned_res[0]/(realigned_res[0]+realigned_res[3])),file=current_file)
            print("illumina 30x\tTP:{}\tTN:{}\tFP:{}\tFN:{}".format(*illumina_res),file=current_file)
            print("precision: {} recall: {}".format(illumina_res[0]/(illumina_res[0]+illumina_res[2]), illumina_res[0]/(illumina_res[0]+illumina_res[3])),file=current_file)


def count_accuracy(realigned_vcf, illumina_vcf, giab_ref_vcf, bed_filter_file, chrom_filter, hg19_fasta, report):

    # files with variants we want to test the accuracy of
    quals = [10,20,30,40,50,60,70,80,90,100,200,300,400]


    # read the variant calls for the reference dataset
    giab_ref = make_vcf_dict(giab_ref_vcf)
    bed_filter = bed_file(bed_filter_file,chrom_filter=chrom_filter)

    # compute confusion matrix for each set of calls
    realigned_precisions = []
    realigned_recalls    = []
    illumina_precisions  = []
    illumina_recalls     = []

    for qual in quals:
        #qual = 30

        print(qual)
        #print("TESTING PACBIO REALIGNMENT:")
        realigned_test = make_vcf_dict(realigned_vcf, min_gq=qual)

        realigned_res, realigned_FP_errors = compare_test_ref(realigned_test,giab_ref,bed_filter, chrom_filter)
        #if ((realigned_res[0] + realigned_res[2]) == 0 or (realigned_res[0]+realigned_res[3]) == 0):
        #    break
        if realigned_res[0] + realigned_res[2] == 0 or realigned_res[0]+realigned_res[3] == 0:
            continue

        p = realigned_res[0] / (realigned_res[0] + realigned_res[2])
        r = realigned_res[0]/(realigned_res[0]+realigned_res[3])
        if p != 0 and r != 0:
            realigned_precisions.append(p)
            realigned_recalls.append(r)

        #print("TESTING STANDARD ILLUMINA:")
        illumina_test = make_vcf_dict(illumina_vcf,min_qual=qual)
        illumina_res, illumina_FP_errors = compare_test_ref(illumina_test,giab_ref,bed_filter, chrom_filter)
        #if ((illumina_res[0] + illumina_res[2]) == 0 or (illumina_res[0]+illumina_res[3]) == 0):
        #    break
        illumina_precisions.append(illumina_res[0] / (illumina_res[0] + illumina_res[2]))
        illumina_recalls.append(illumina_res[0]/(illumina_res[0]+illumina_res[3]))

        #print("")
        #if qual == quals[5]:
        with open(report,'w') as outf:
            for current_file in [outf,sys.stdout]:
                print("pacbio realignment\tTP:{}\tTN:{}\tFP:{}\tFN:{}".format(*realigned_res),file=current_file)
                print("precision: {} recall: {}".format(realigned_res[0]/(realigned_res[0]+realigned_res[2]), realigned_res[0]/(realigned_res[0]+realigned_res[3])),file=current_file)
                print("illumina 30x\tTP:{}\tTN:{}\tFP:{}\tFN:{}".format(*illumina_res),file=current_file)
                print("precision: {} recall: {}".format(illumina_res[0]/(illumina_res[0]+illumina_res[2]), illumina_res[0]/(illumina_res[0]+illumina_res[3])),file=current_file)

    prec_recalls = realigned_precisions, realigned_recalls, illumina_precisions, illumina_recalls
    pickle.dump(prec_recalls,open('prec_recalls.p','wb'))
    realigned_precisions, realigned_recalls, illumina_precisions, illumina_recalls = pickle.load(
        open('prec_recalls.p', 'rb'))
    plt.figure()
    plt.plot(realigned_recalls, realigned_precisions, color='r', label='Realigned PacBio')
    plt.plot(illumina_recalls, illumina_precisions, color='b', label='Illumina 30x')
    plt.xlim((0.5,1.005))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Realigned PacBio vs Illumina Genotyping:\nAll GIAB Confident Regions')
    plt.legend()
    plt.savefig('giab_confident_all.png')
