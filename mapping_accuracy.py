
import pysam
import re
from collections import defaultdict
from matplotlib import pyplot as plt

# bamfile should be filtered to only be segdups
# bedfile contains the segmental duplications
# chrom is an optional chromosome filter
def mapping_accuracy_and_completeness_segdups(bamfile, bedfile, outfile, chrom_filter=None, min_mapq=30, delta=5000):

    r = re.compile('startpos=(\d+)')

    with pysam.AlignmentFile(bamfile, "rb") as bam, open(bedfile, 'r') as bed, open(outfile,'w') as outf:
        for line in bed:

            # parse the bedfile with segmental duplication regions
            el = line.strip().split()
            chrom = el[0]
            start = int(el[1]) # 0-based
            stop = int(el[2])  # 0-based, partially open

            if chrom_filter != None and chrom != chrom_filter:
                continue

            # count the mean coverage in this region
            cov_total = 0
            pos10 = 0
            pos20 = 0
            pos30 = 0
            pos40 = 0
            pos50 = 0
            pos60 = 0

            for pileupcolumn in bam.pileup(chrom, start, stop, truncate=True):
                n = pileupcolumn.nsegments
                if n >= 10:
                    pos10 += 1
                if n >= 20:
                    pos20 += 1
                if n >= 30:
                    pos30 += 1
                if n >= 40:
                    pos40 += 1
                if n >= 50:
                    pos50 += 1
                if n >= 60:
                    pos60 += 1

                cov_total += n

            region_size = (stop - start)
            assert(region_size > 0)
            
            map_coverage = cov_total / region_size
            cov10_frac = pos10 / region_size
            cov20_frac = pos20 / region_size
            cov30_frac = pos30 / region_size
            cov40_frac = pos40 / region_size
            cov50_frac = pos50 / region_size
            cov60_frac = pos60 / region_size

            # calculate the mapping accuracy in this segmental dup

            correct = 0
            incorrect = 0
            #incorrect_pos_lst = []

            for record in bam.fetch(chrom, start, stop):

                if record.mapq < min_mapq:
                    continue

                m = r.search(record.qname)
                correct_pos = int(m.group(1))
                if abs(record.pos - correct_pos) < delta:
                    correct += 1
                else:
                    incorrect += 1
                    #incorrect_pos_lst.append(record.pos)

            if  (correct + incorrect) > 0:
                mapping_accuracy = correct / (correct + incorrect)
            else:
                mapping_accuracy = 0

            # print line
            line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom,start,stop,map_coverage,mapping_accuracy,
            cov10_frac,cov20_frac,cov30_frac,cov40_frac,cov50_frac,cov60_frac)
            print(line, file = outf)

def mapping_accuracy(bamfile, plot_name, min_mapq=30, delta=5000):

    correct = 0
    incorrect = 0
    incorrect_pos_lst = []

    r = re.compile('startpos=(\d+)')
    with pysam.AlignmentFile(bamfile, "rb") as bf:

        for record in bf:

            if record.mapq < min_mapq:
                continue

            m = r.search(record.qname)
            correct_pos = int(m.group(1))
            if abs(record.pos - correct_pos) < delta:
                correct += 1
            else:
                incorrect += 1
                incorrect_pos_lst.append(record.pos)

    accuracy = correct / (correct + incorrect)
    print("mapping accuracy: {}".format(accuracy))

    ################################################
    incorrect_pos_lst.sort()
    prev_pos = None
    binsize = 100000
    bin_dict = defaultdict(int)

    for pos in incorrect_pos_lst:

        binned_pos = int(pos / binsize)
        bin_dict[binned_pos] += 1

    counts = []
    #for x in range(0,max(bin_dict.keys())):
    for x in bin_dict.keys():
        counts.append(bin_dict[x])

    fig = plt.figure()
    ax = plt.subplot(111)
    plt.hist(counts, alpha=0.6)
    #ax.set_yscale('log')
    # axis formatting

    # hiding axis ticks
    plt.tick_params(axis="both", which="both", bottom=False, top=False,
            labelbottom=True, left=False, right=False, labelleft=True)

    # adding horizontal grid lines
    ax.yaxis.grid(True,linestyle='--',color='grey',alpha=0.5)

    # remove axis spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.title("Distribution of Mismapped Reads at Loci with >=1 Mismapping")
    plt.xlabel("Number of Mismapped Reads")
    plt.ylabel("Number of Loci")
    plt.savefig(plot_name)
    ################################################

    return accuracy

#mapping_accuracy('data/simulation/aligned_reads/pacbio/pacbio.bwa.1.20x.segdup.bam','plt.png')
