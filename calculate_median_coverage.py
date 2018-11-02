
import pysam
import statistics
import time

CHROM_FILTER = set(['chr{}'.format(x) for x in range(1,23)]+['{}'.format(x) for x in range(1,23)])


def calculate_median_coverage(bam_file, random_positions_bed, chrom_filter=CHROM_FILTER, min_mapq=30, flag_filter=3844, add_chr=False):

    t1 = time.time()

    depth_lst = []

    with pysam.AlignmentFile(bam_file, "rb") as bam, open(random_positions_bed,'r') as bed:

        for line in bed:

            # parse the bedfile with segmental duplication regions
            el = line.strip().split()
            chrom = el[0]
            start = int(el[1]) # 0-based
            stop = int(el[2])  # 0-based, partially open
            if add_chr:
                chrom = 'chr' + chrom

            if chrom not in chrom_filter:
                continue

            assert(stop == start + 1)

            # should only check one position
            for pileup_column in bam.pileup(chrom, start, stop, truncate=True,
                                           flag_filter=flag_filter, min_mapping_quality=min_mapq, min_base_quality=0):
                assert(pileup_column.reference_pos == start)
                dp = 0# 0
                for pileup_read in pileup_column.pileups:
                    dp += 1
                #assert(dp == pileup_column.nsegments)

                depth_lst.append(dp)

    t2 = time.time()
    print("time: {:2}".format(t2-t1))

    return statistics.median(depth_lst)
