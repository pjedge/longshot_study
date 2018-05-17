
import pysam
from collections import defaultdict

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
                     line_tabix_bed_file, sine_tabix_bed_file, ref_fa):

    counts = defaultdict(int)
    total = 0

    with pysam.VariantFile(calls_vcfgz) as calls, \
        pysam.VariantFile(ground_truth_vcfgz) as ground_truth, \
        pysam.TabixFile(str_tabix_bed_file) as str_bed, \
        pysam.TabixFile(line_tabix_bed_file) as line_bed, \
        pysam.TabixFile(sine_tabix_bed_file) as sine_bed, \
        pysam.FastaFile(ref_fa) as ref:

        for call in calls:

            total += 1
            assert(call.chrom == chrom_name)
            # does the variant border on a homopolymer of length 5?
            has_homopol5 = has_homopolymer(ref, call.chrom, call.pos, 5, 5)

            # does the variant occur within 30 bp of a true indel?
            near_indel = False
            indel_pad = 30
            for rec in ground_truth.fetch(contig=call.chrom,start=call.pos-indel_pad-1,stop=call.pos+indel_pad):
                is_indel = (len(rec.samples[0].alleles[0]) != len(rec.ref) or
                            len(rec.samples[0].alleles[1]) != len(rec.ref))

                if is_indel:
                    near_indel = True
                    break

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

            observed_categories = 0
            # iterate appropriate counts
            if has_homopol5:
                counts['Inside homopolymer length >= 5:'] += 1
                observed_categories += 1
            if near_indel:
                counts['Near indel:'] += 1
                observed_categories += 1
            if in_STR:
                counts['Inside STR:'] += 1
                observed_categories += 1
            if in_LINE:
                counts['Inside LINE:'] += 1
                observed_categories += 1
            if in_SINE:
                counts['Inside SINE:'] += 1
                observed_categories += 1

            if observed_categories >= 2:
                counts['Multiple categories:'] += 1

    for category, count in counts.items():
        if category == 'Multiple categories:':
            continue

        # pad the category string
        assert(len(category) <= 40)
        category += ''.join([' ']*(40-len(category)))

        print("{} {:.3f}".format(category, count / total))

    print('\nNote that these categories are not mutually exclusive and there may be overlap!')
    print('{} {:.3f}'.format('Multiple categories:', counts['Multiple categories:'] / total))
