
import pysam
from collections import defaultdict, namedtuple
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from random import randrange

category_data = namedtuple('category_data',['near_indel', 'in_homopol5', 'in_STR', 'in_LINE', 'in_SINE'])
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

def generate_random_calls(ground_truth_bed_file, chrom, N):

    var_pos_lst = []
    with pysam.TabixFile(ground_truth_bed_file) as bed:
        end = 0
        for record in bed.fetch(reference=chrom,parser=pysam.asBed()):
            end = record.end

        while len(var_pos_lst) < N:
            randpos = randrange(0, end+1)
            # does the variant occur within 5 bases of an STR?
            in_confident_region = False
            for row in bed.fetch(chrom, randpos-1, randpos, parser=pysam.asBed()):
                in_confident_region = True

            if in_confident_region:
                var_pos_lst.append((chrom, randpos))

    return var_pos_lst

def get_var_pos_lst(calls_vcfgz):

    var_pos_lst = []

    with pysam.VariantFile(calls_vcfgz) as calls:

        for record in calls:
            var_pos_lst.append((record.chrom, record.pos))

    return var_pos_lst

def analyze_variants(chrom_name, fp_calls_vcfgz, fn_calls_vcfgz, ground_truth_vcfgz, ground_truth_bed_file, str_tabix_bed_file,
                     line_tabix_bed_file, sine_tabix_bed_file, ref_fa, output_file):

    def count_variant_categories(var_pos_lst):

        counts = defaultdict(int)

        total = 0
        near_indel_ct = 0
        in_homopol5_ct = 0
        in_STR_ct = 0
        in_LINE_ct = 0
        in_SINE_ct = 0

        with pysam.VariantFile(ground_truth_vcfgz) as ground_truth, \
            pysam.TabixFile(str_tabix_bed_file) as str_bed, \
            pysam.TabixFile(line_tabix_bed_file) as line_bed, \
            pysam.TabixFile(sine_tabix_bed_file) as sine_bed, \
            pysam.FastaFile(ref_fa) as ref:

            for (chrom, pos) in var_pos_lst:

                total += 1
                assert(chrom == chrom_name)

                # does the variant occur within 30 bp of a true indel?
                near_indel = False
                indel_pad = 30
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
                str_pad = 5
                for row in str_bed.fetch(chrom, pos-1-str_pad, pos+str_pad, parser=pysam.asBed()):
                    in_STR = True

                # does the variant occur within 5 bases of a LINE?
                in_LINE = False
                line_pad = 5
                for row in line_bed.fetch(chrom, pos-1-line_pad, pos+line_pad, parser=pysam.asBed()):
                    in_LINE = True

                # does the variant occur within 5 bases of a SINE?
                in_SINE = False
                sine_pad = 5
                for row in sine_bed.fetch(chrom, pos-1-sine_pad, pos+sine_pad, parser=pysam.asBed()):
                    in_SINE = True

                near_indel_ct += near_indel
                in_homopol5_ct += in_homopol5
                in_STR_ct += in_STR
                in_LINE_ct += in_LINE
                in_SINE_ct += in_SINE

                counts[category_data(near_indel=near_indel, in_homopol5=in_homopol5,
                       in_STR=in_STR, in_LINE=in_LINE, in_SINE=in_SINE)] += 1

        assert(total == sum(counts.values()))

        #bit_table = []
        #row_labels = []
        result_fracs = [near_indel_ct/total,
                             in_homopol5_ct/total,
                             in_STR_ct/total,
                             in_LINE_ct/total,
                             in_SINE_ct/total]

        print("Fraction of variants in categories (categories may overlap):")
        print("Near true indel (within 30 bp): {:.3f}".format(near_indel_ct/total))
        print("In homopolymer (len >= 5):      {:.3f}".format(in_homopol5_ct/total))
        print("In STR (+- 5 bp):               {:.3f}".format(in_STR_ct/total))
        print("In LINE (+- 5 bp):              {:.3f}".format(in_LINE_ct/total))
        print("In SINE (+- 5 bp):              {:.3f}".format(in_SINE_ct/total))
        print("")
        print("")
        print("Fractions for overlapped categories:")
        print("Near indel\tIn homopolymer\tIn STR\tIn LINE\tIn SINE\tFraction of Variants")
        sorted_counts = sorted(counts.items(),key=lambda x: x[0])

        for cats, count in sorted_counts:

            print("{}\t{}\t{}\t{}\t{}\t{:.3f}".format(int(cats.near_indel), int(cats.in_homopol5),
            int(cats.in_STR), int(cats.in_LINE), int(cats.in_SINE), count / total))
            #bit_table.append((int(cats.near_indel), int(cats.in_homopol5),
            #int(cats.in_STR), int(cats.in_LINE), int(cats.in_SINE)))
            #row_labels.append("{:.3f}".format(count / total))

        #return bit_table, row_labels, col_labels_bottom
        return result_fracs

    random_var_pos = generate_random_calls(ground_truth_bed_file, chrom_name, 100000)

    print("Analyzing False Positives...\n")
    fp_fracs = count_variant_categories(get_var_pos_lst(fp_calls_vcfgz))
    print("Analyzing False Negatives...\n")
    fn_fracs = count_variant_categories(get_var_pos_lst(fn_calls_vcfgz))
    print("Analyzing Random Positions...\n")
    random_fracs = count_variant_categories(random_var_pos)
    print(fp_fracs)
    print(fn_fracs)
    print(random_fracs)

    s = '''
\\begin{{table}}[htbp]
\\centering
\\begin{{tabular}}{{lrrrrr}}
\\hline
Genome      & False          & False          & Random    & FP/random & FN/random  \\\\
            & Positives (FP) & Negatives (FN) & Positions &           &            \\\\
\\hline
Near indel     & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} \\\\
In homopolymer & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} \\\\
In STR         & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} \\\\
In LINE        & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} \\\\
In SINE        & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} \\\\
\\hline
\\end{{tabular}}
\\caption{{{{\\bf Fractions of False Positive (FN) and False Negative (FN) variant calls that
coincide with genomic locations or features. For comparison, random positions from
the GIAB confident regions were selected and subjected to the same analysis. The
last two columns show the fold difference compared to the random positions.}}}}
\\label{{tab:stats}}
\\end{{table}}
'''.format(fp_fracs[0], fn_fracs[0], random_fracs[0], fp_fracs[0]/random_fracs[0], fn_fracs[0]/random_fracs[0],
           fp_fracs[1], fn_fracs[1], random_fracs[1], fp_fracs[1]/random_fracs[1], fn_fracs[1]/random_fracs[1],
           fp_fracs[2], fn_fracs[2], random_fracs[2], fp_fracs[2]/random_fracs[2], fn_fracs[2]/random_fracs[2],
           fp_fracs[3], fn_fracs[3], random_fracs[3], fp_fracs[3]/random_fracs[3], fn_fracs[3]/random_fracs[3],
           fp_fracs[4], fn_fracs[4], random_fracs[4], fp_fracs[4]/random_fracs[4], fn_fracs[4]/random_fracs[4])

    with open(output_file,'w') as outf:
        print(s, file=outf)

    col_labels = ["Near indel", "In homopolymer", "In STR", "In LINE", "In SINE"]

    '''
    # credit to this SO thread: https://stackoverflow.com/questions/33283601/manually-defined-axis-labels-for-matplotlib-imshow

    ###########################################################################
    # plot false positive table
    ###########################################################################

    fig = plt.figure()
    ax1 = SubplotHost(fig, 1,2,1)
    fig.add_subplot(ax1)

    # set bottom ticks
    ax1.set_xticks(list(range(0,len(col_labels))))
    ax1.set_xticklabels(fp_col_labels_bottom,rotation='vertical')
    #ax1.set_xlabel("Total Fraction\n in Category")
    ax1.tick_params(length=0)
    ax1.imshow(fp_bit_table, cmap='Greys', interpolation='nearest')

    ax2 = ax1.twin()
    # set top ticks
    ax2.set_xticks(list(range(0,len(col_labels))))
    ax2.set_xticklabels(col_labels,rotation='vertical')
    ax2.xaxis.tick_top()
    ax2.tick_params(length=0)
    ax2.set_yticks([])

    # set row ticks
    ax1.set_yticks(list(range(0,len(fp_row_labels))))
    ax1.set_yticklabels(fp_row_labels)
    plt.ylabel("Fraction of variants\n in intersection of categories")

    plt.title("False Positive Variants",y=1.5)
    #plt.subplots_adjust(top=0.8, bottom=0.15)

    ###########################################################################
    # plot false negative table
    ###########################################################################

    ax1 = SubplotHost(fig, 1,2,2)

    fig.add_subplot(ax1)

    # set bottom ticks
    ax1.set_xticks(list(range(0,len(col_labels))))
    ax1.set_xticklabels(fn_col_labels_bottom,rotation='vertical')
    #ax1.set_xlabel("Total Fraction\n in Category")
    ax1.tick_params(length=0)
    ax1.imshow(fn_bit_table, cmap='Greys', interpolation='nearest')

    ax2 = ax1.twin()
    # set top ticks
    ax2.set_xticks(list(range(0,len(col_labels))))
    ax2.set_xticklabels(col_labels,rotation='vertical')
    ax2.xaxis.tick_top()
    ax2.tick_params(length=0)
    ax2.set_yticks([])

    # set row ticks
    ax1.set_yticks(list(range(0,len(fn_row_labels))))
    ax1.set_yticklabels(fn_row_labels)
    #plt.ylabel("Fraction of variants\n in intersection of categories")
    plt.title("False Negative Variants",y=1.5)
    fig.text(0.5, -0.03, "Total Fraction of Variants in Category", ha='center')

    plt.subplots_adjust(top=0.7, bottom=0.15, wspace=0.2)
    plt.savefig(output_figure)
    '''
