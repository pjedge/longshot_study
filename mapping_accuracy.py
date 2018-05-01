
import pysam
import re
from collections import defaultdict
from matplotlib import pyplot as plt

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
