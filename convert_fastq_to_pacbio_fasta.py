import pysam

def convert_fastq_to_pacbio_fasta(input_fastq, output_fasta):
    with pysam.FastxFile(input_fastq) as inf, open(output_fasta,'w') as outf:
        for i,record in enumerate(inf):
            print(">{}/{0}/0_{}".format(record.name,i,len(record.seq)),file=outf)
            print(record.seq, file=outf)
