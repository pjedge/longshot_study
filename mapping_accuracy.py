
import pysam
import re

def mapping_accuracy(bamfile, min_mapq=30, delta=5000):

    correct = 0
    incorrect = 0
    r = re.compile('startpos=(\d+)')
    with pysam.AlignmentFile(bamfile, "rb") as bf:

        for record in bf:

            if record.mapq < min_mapq:
                continue

            m = r.search(record.qname)
            if abs(record.pos - int(m.group(1))) < delta:
                correct += 1
            else:
                incorrect += 1

    accuracy = correct / (correct + incorrect)
    print("mapping accuracy: {}".format(accuracy))

    return accuracy

mapping_accuracy('data/simulation/aligned_reads/pacbio/pacbio.bwa.1.20x.bam')
