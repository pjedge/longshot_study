
import pysam
from collections import defaultdict, namedtuple
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from random import randrange

category_data = namedtuple('category_data',['misgenotyped','near_indel', 'in_homopol5', 'in_STR', 'in_LINE', 'in_SINE'])
# INPUT
# fasta_object: a pysam indexed fasta object
# chrom, pos: chromosome name and 1-indexed position to analyze
# pad: number of bases to the left and right of (chrom, pos) to consider in the window
# run_length: length of homopolymer run to count as a homopolymer

# OUTPUT
# boolean value, whether reference has a homopolymer at this position based on the criteria
def has_homopolymer(fasta_object, chrom, pos, pad=5, run_length=3):

    window = str.upper(fasta_object.fetch(chrom, pos-pad-1, pos+pad))
    curr_letter = window[0]
    count = 1

    for letter in window[1:]:
        if letter == curr_letter:
            count += 1
        else:
            count = 1
            curr_letter = letter

        if count >= run_length:
            return True

    return False

def generate_random_calls(chrom, chrlen, N, outfile):

    pos_lst = []

    for i in range(N):
        pos_lst.append(randrange(0, chrlen))

    pos_lst.sort()

    with open(outfile, 'w') as outf:

        header = '''##fileformat=VCFv4.2
##source=analyze_variants.py
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRANDOM'''
        print(header, file=outf)
        for pos in pos_lst:
            print("{}\t{}\t.\tN\tN\t100\tPASS\t.\tGT:GQ\t1/1:100".format(chrom, pos),file=outf)

def get_var_pos_lst(calls_vcfgz, gq_cutoff):

    var_pos_lst = []

    with pysam.VariantFile(calls_vcfgz) as calls:

        for record in calls:
            if record.samples[0]['GQ'] < gq_cutoff:
                continue

            # currently only designed for variants that are simple SNVs
            assert(len(record.alleles[0]) == 1)

            var_pos_lst.append((record.chrom, record.pos, record.alleles[0],
                                record.samples[0].alleles))

    return var_pos_lst

def analyze_variants(chrom_name, pacbio_calls_vcfgz, fp_calls_vcfgz, fn_calls_vcfgz, ground_truth_vcfgz, ground_truth_bed_file, random_positions_vcfgz,
                     str_tabix_bed_file, line_tabix_bed_file, sine_tabix_bed_file, ref_fa, gq_cutoff, output_file):

    def count_variant_categories(var_pos_lst, mode):

        assert(mode in ['fp','fn','rand'])
        counts = defaultdict(int)

        total = 0
        misgenotyped_ct = 0
        near_indel_ct = 0
        in_homopol5_ct = 0
        in_homopol5_not_near_indel_ct = 0
        in_STR_ct = 0
        in_LINE_ct = 0
        in_SINE_ct = 0

        with pysam.VariantFile(pacbio_calls_vcfgz) as pacbio_calls, \
            pysam.VariantFile(ground_truth_vcfgz) as ground_truth, \
            pysam.TabixFile(str_tabix_bed_file) as str_bed, \
            pysam.TabixFile(line_tabix_bed_file) as line_bed, \
            pysam.TabixFile(sine_tabix_bed_file) as sine_bed, \
            pysam.FastaFile(ref_fa) as ref:

            # r, v are strings representing ref allele base, alt base
            # gt is tuple of ints representing genotype
            for ix, (chrom, pos, ref_base, alleles) in enumerate(var_pos_lst):

                total += 1
                assert(chrom == chrom_name)

                # was the variant found, but misgenotyped?
                misgenotyped = False
                if mode in ['fp','fn']:
                    # if we're analyzing FPs then in the FP file we had a list of
                    # pacbio calls. we want to compare to the ground truth calls
                    # if we're analyzing FNs then in the FN file we had a list of
                    # ground truth calls. we want to compare to the pacbio calls.
                    calls = ground_truth if mode == 'fp' else pacbio_calls
                    for rec in calls.fetch(contig=chrom,start=pos-1,stop=pos):
                        if rec.pos != pos:
                            continue

                        SNV = True
                        for a in rec.alleles:
                            if len(a) != 1:
                                SNV = False

                        if not SNV:
                            continue

                        assert(ref_base == rec.ref)

                        if not set(rec.samples[0].alleles) == set(alleles):
                            misgenotyped = True
                            break

                if misgenotyped:
                    print("{} {} was misgenotyped".format(chrom, pos))

                # does the variant occur within 30 bp of a true indel?
                near_indel = False
                indel_pad = 10
                for rec in ground_truth.fetch(contig=chrom,start=pos-indel_pad-1,stop=pos+indel_pad):
                    is_indel = (len(rec.samples[0].alleles[0]) != len(rec.ref) or
                                len(rec.samples[0].alleles[1]) != len(rec.ref))

                    if is_indel:
                        near_indel = True
                        break

                # does the variant border on a homopolymer of length 5?
                in_homopol5 = has_homopolymer(ref, chrom, pos, 5, 5)

                # does the variant occur within 5 bases of an STR?
                in_STR = False
                str_pad = 0
                for row in str_bed.fetch(chrom, pos-1-str_pad, pos+str_pad, parser=pysam.asBed()):
                    in_STR = True

                # does the variant occur within 5 bases of a LINE?
                in_LINE = False
                line_pad = 0
                for row in line_bed.fetch(chrom, pos-1-line_pad, pos+line_pad, parser=pysam.asBed()):
                    in_LINE = True

                # does the variant occur within 5 bases of a SINE?
                in_SINE = False
                sine_pad = 0
                for row in sine_bed.fetch(chrom, pos-1-sine_pad, pos+sine_pad, parser=pysam.asBed()):
                    in_SINE = True

                misgenotyped_ct += misgenotyped
                near_indel_ct += near_indel
                in_homopol5_ct += in_homopol5
                in_homopol5_not_near_indel_ct += in_homopol5 and not near_indel
                in_STR_ct += in_STR
                in_LINE_ct += in_LINE
                in_SINE_ct += in_SINE

                counts[category_data(misgenotyped=misgenotyped, near_indel=near_indel, in_homopol5=in_homopol5,
                                     in_STR=in_STR, in_LINE=in_LINE, in_SINE=in_SINE)] += 1

        assert(total == sum(counts.values()))

        #bit_table = []
        #row_labels = []
        result_fracs = [misgenotyped_ct/total,
                        near_indel_ct/total,
                        in_homopol5_ct/total,
                        in_homopol5_not_near_indel_ct/total,
                        in_STR_ct/total,
                        in_LINE_ct/total,
                        in_SINE_ct/total]


        print("Counts of variants in categories (categories may overlap):")
        print("Misgenotyped:                   {:.3f}".format(misgenotyped_ct))
        print("Near true indel (within 10 bp): {:.3f}".format(near_indel_ct))
        print("In homopolymer (len >= 5):      {:.3f}".format(in_homopol5_ct))
        print("In homopolymer, not near indel: {:.3f}".format(in_homopol5_not_near_indel_ct))
        print("In STR:               {:.3f}".format(in_STR_ct))
        print("In LINE:              {:.3f}".format(in_LINE_ct))
        print("In SINE:              {:.3f}".format(in_SINE_ct))
        print("")
        print("")
        print("Fractions for overlapped categories:")
        print("Near indel\tIn homopolymer\tIn STR\tIn LINE\tIn SINE\tFraction of Variants")

        print("Fraction of variants in categories (categories may overlap):")
        print("Misgenotyped:                   {:.3f}".format(misgenotyped_ct/total))
        print("Near true indel (within 10 bp): {:.3f}".format(near_indel_ct/total))
        print("In homopolymer (len >= 5):      {:.3f}".format(in_homopol5_ct/total))
        print("In homopolymer, not near indel: {:.3f}".format(in_homopol5_not_near_indel_ct/total))
        print("In STR:               {:.3f}".format(in_STR_ct/total))
        print("In LINE:              {:.3f}".format(in_LINE_ct/total))
        print("In SINE:              {:.3f}".format(in_SINE_ct/total))
        print("")
        print("")
        print("Fractions for overlapped categories:")
        print("Misgenotyped\tNear indel\tIn homopolymer\tIn STR\tIn LINE\tIn SINE\tFraction of Variants")
        sorted_counts = sorted(counts.items(),key=lambda x: x[0])

        for cats, count in sorted_counts:

            print("{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}".format(int(cats.misgenotyped),int(cats.near_indel), int(cats.in_homopol5),
            int(cats.in_STR), int(cats.in_LINE), int(cats.in_SINE), count / total))

        return result_fracs

    #N = 100000
    #random_var_pos = generate_random_calls(ground_truth_bed_file, chrom_name, N)
    print("Analyzing False Positives...\n")
    fp_fracs = count_variant_categories(get_var_pos_lst(fp_calls_vcfgz, gq_cutoff), 'fp')
    print("Analyzing False Negatives...\n")
    fn_fracs = count_variant_categories(get_var_pos_lst(fn_calls_vcfgz, gq_cutoff), 'fn')
    print("Analyzing Random Positions...\n")
    random_fracs = count_variant_categories(get_var_pos_lst(random_positions_vcfgz, gq_cutoff), 'rand')
    print(fp_fracs)
    print(fn_fracs)
    print(random_fracs)

    #for i in range(0,len(random_fracs)):
    #    if random_fracs[i] == 0:
    #        random_fracs[i] = (1/N)**2

    s = '''
\\begin{{table}}[htbp]
\\centering
\\begin{{tabular}}{{lllll}}
\\hline
Genome      & False          & FP         & False          & FN         \\\\
            & Positives (FP) & Enrichment & Negatives (FN) & Enrichment \\\\
\\hline
Misgenotyped SNV & {:.3f} & - & {:.3f} & - \\\\
Near Indel     & {:.3f} & {:.2f} & {:.3f} & {:.2f} \\\\
In homopolymer & {:.3f} & {:.2f} & {:.3f} & {:.2f} \\\\
In homopolymer but & {:.3f} & {:.2f} & {:.3f} & {:.2f} \\\\
\\ \\ \\ \\ \\ not near indel &      &        &        &        \\\\
In STR         & {:.3f} & {:.2f} & {:.3f} & {:.2f} \\\\
In LINE        & {:.3f} & {:.2f} & {:.3f} & {:.2f} \\\\
In SINE        & {:.3f} & {:.2f} & {:.3f} & {:.2f} \\\\
\\hline
\\end{{tabular}}
\\caption{{{{\\bf Fractions of False Positive (FN) and False Negative (FN) variant calls that
were misgenotyped or coincide with genomic features. For comparison, random positions from
the GIAB confident regions were selected and subjected to the same analysis. The
last two columns show the fold enrichment compared to the random positions.}}}}
\\label{{tab:stats}}
\\end{{table}}
'''.format(fp_fracs[0], fn_fracs[0],
           fp_fracs[1], fp_fracs[1]/random_fracs[1], fn_fracs[1], fn_fracs[1]/random_fracs[1],
           fp_fracs[2], fp_fracs[2]/random_fracs[2], fn_fracs[2], fn_fracs[2]/random_fracs[2],
           fp_fracs[3], fp_fracs[3]/random_fracs[3], fn_fracs[3], fn_fracs[3]/random_fracs[3],
           fp_fracs[4], fp_fracs[4]/random_fracs[4], fn_fracs[4], fn_fracs[4]/random_fracs[4],
           fp_fracs[5], fp_fracs[5]/random_fracs[5], fn_fracs[5], fn_fracs[5]/random_fracs[5],
           fp_fracs[6], fp_fracs[6]/random_fracs[6], fn_fracs[6], fn_fracs[6]/random_fracs[6])

    with open(output_file,'w') as outf:
        print(s, file=outf)

def count_fp_near_true_indel(fp_calls_vcfgz, ground_truth_vcfgz, gq_cutoff):

    near_indel_ct = 0
    indel_pad = 10
    var_pos_lst = get_var_pos_lst(fp_calls_vcfgz, gq_cutoff)

    with pysam.VariantFile(ground_truth_vcfgz) as ground_truth:

        for ix, (chrom, pos, ref_base, alleles) in enumerate(var_pos_lst):
            # does the variant occur within 30 bp of a true indel?
            near_indel = False
            for rec in ground_truth.fetch(contig=chrom,start=pos-indel_pad-1,stop=pos+indel_pad):
                is_indel = (len(rec.samples[0].alleles[0]) != len(rec.ref) or
                            len(rec.samples[0].alleles[1]) != len(rec.ref))

                if is_indel:
                    near_indel = True
                    break

            near_indel_ct += near_indel

    return near_indel_ct
