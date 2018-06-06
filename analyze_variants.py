
import pysam
from collections import defaultdict, namedtuple

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

def analyze_variants(chrom_name, calls_vcfgz, ground_truth_vcfgz, str_tabix_bed_file,
                     line_tabix_bed_file, sine_tabix_bed_file, ref_fa, output_file):

    counts = defaultdict(int)

    total = 0
    near_indel_ct = 0
    in_homopol5_ct = 0
    in_STR_ct = 0
    in_LINE_ct = 0
    in_SINE_ct = 0

    with pysam.VariantFile(calls_vcfgz) as calls, \
        pysam.VariantFile(ground_truth_vcfgz) as ground_truth, \
        pysam.TabixFile(str_tabix_bed_file) as str_bed, \
        pysam.TabixFile(line_tabix_bed_file) as line_bed, \
        pysam.TabixFile(sine_tabix_bed_file) as sine_bed, \
        pysam.FastaFile(ref_fa) as ref:

        for call in calls:

            total += 1
            assert(call.chrom == chrom_name)

            # does the variant occur within 30 bp of a true indel?
            near_indel = False
            indel_pad = 30
            for rec in ground_truth.fetch(contig=call.chrom,start=call.pos-indel_pad-1,stop=call.pos+indel_pad):
                is_indel = (len(rec.samples[0].alleles[0]) != len(rec.ref) or
                            len(rec.samples[0].alleles[1]) != len(rec.ref))

                if is_indel:
                    near_indel = True
                    break

            # does the variant border on a homopolymer of length 5?
            in_homopol5 = has_homopolymer(ref, call.chrom, call.pos, 5, 5)

            # does the variant occur within 5 bases of an STR?
            in_STR = False
            str_pad = 5
            for row in str_bed.fetch(call.chrom, call.pos-1-str_pad, call.pos+str_pad, parser=pysam.asBed()):
                in_STR = True

            # does the variant occur within 5 bases of a LINE?
            in_LINE = False
            line_pad = 5
            for row in line_bed.fetch(call.chrom, call.pos-1-line_pad, call.pos+line_pad, parser=pysam.asBed()):
                in_LINE = True

            # does the variant occur within 5 bases of a SINE?
            in_SINE = False
            sine_pad = 5
            for row in sine_bed.fetch(call.chrom, call.pos-1-sine_pad, call.pos+sine_pad, parser=pysam.asBed()):
                in_SINE = True

            near_indel_ct += near_indel
            in_homopol5_ct += in_homopol5
            in_STR_ct += in_STR
            in_LINE_ct += in_LINE
            in_SINE_ct += in_SINE

            counts[category_data(near_indel=near_indel, in_homopol5=in_homopol5,
                   in_STR=in_STR, in_LINE=in_LINE, in_SINE=in_SINE)] += 1

    assert(total == sum(counts.values()))

    with open(output_file, 'w') as outf:

        print("Fraction of variants in categories (categories may overlap):", file=outf)
        print("Near true indel (within 30 bp): {:.3f}".format(near_indel_ct/total), file=outf)
        print("In homopolymer (len >= 5):      {:.3f}".format(in_homopol5_ct/total), file=outf)
        print("In STR (+- 5 bp):               {:.3f}".format(in_STR_ct/total), file=outf)
        print("In LINE (+- 5 bp):              {:.3f}".format(in_LINE_ct/total), file=outf)
        print("In SINE (+- 5 bp):              {:.3f}".format(in_SINE_ct/total), file=outf)
        print("", file=outf)
        print("", file=outf)
        print("Fractions for overlapped categories:", file=outf)
        print("Near indel\tIn homopolymer\tIn STR\tIn LINE\tIn SINE\tFraction of Variants", file=outf)
        sorted_counts = sorted(counts.items(),key=lambda x: x[0])

        for cats, count in sorted_counts:

            print("{}\t{}\t{}\t{}\t{}\t{:.3f}".format(int(cats.near_indel), int(cats.in_homopol5),
            int(cats.in_STR), int(cats.in_LINE), int(cats.in_SINE), count / total),file=outf)
