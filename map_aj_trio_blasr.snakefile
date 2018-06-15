
from collections import defaultdict

FASTQUTILS = '/home/pedge/installed/ngsutils/bin/fastqutils'
TWOBITTOFASTA = 'twoBitToFa' # can be downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'
SAMTOOLS       = 'samtools' # v1.3
FASTQ_DUMP     = 'fastq-dump' # v2.5.2
REAPER         = '../target/release/reaper' # v0.1
RTGTOOLS       = 'rtg' #'/home/pedge/installed/rtg-tools-3.8.4/rtg' # v3.8.4, https://www.realtimegenomics.com/products/rtg-tools
BGZIP = 'bgzip'
TABIX = 'tabix'
FREEBAYES      = '/home/pedge/git/freebayes/bin/freebayes'
SIMLORD = '/home/pedge/installed/opt/python/bin/simlord'
DWGSIM = '/home/pedge/git/DWGSIM/dwgsim'
BWA = '/home/pedge/installed/bwa'
BLASR = 'blasr'
BAX2BAM = 'bax2bam'
SAWRITER = 'sawriter'
NGMLR = 'ngmlr'
MINIMAP2 = 'minimap2'
BCFTOOLS = '/opt/biotools/bcftools/bin/bcftools'
PYFAIDX = '/home/pedge/installed/opt/python/bin/faidx'
chroms = ['{}'.format(i) for i in range(1,23)] + ['X']
BEDTOOLS = 'bedtools' # v 2.27

rule all:
    input: expand('study/data/{individual}/aligned_reads/pacbio/pacbio.blasr.all.giab_full_coverage.bam',individual=['NA24143','NA24149'])

# file with 3 columns:
# 1. URL of compressed pacbio hdf5 data
# 2. md5 sum
# 3. HG00X id
GIAB_AJTRIO_PACBIO_HDF5_INDEX = 'sequence.index.AJtrio_PacBio_MtSinai_NIST_hdf5_08072015.txt'

# map the GIAB individuals' names to the number of entries that they should have in the file
num_hdf5_archives = {'NA24385': 292, 'NA24149': 139, 'NA24143': 132}
hg2na = {'HG002':'NA24385', 'HG003':'NA24149', 'HG004':'NA24143'}

# build a dictionary keyed on the individual's NAXXXXX ID
# goes to a list of tuples [(data_url, md5_sum)]
hdf5_file_list = defaultdict(list)
with open(GIAB_AJTRIO_PACBIO_HDF5_INDEX,'r') as infile:
    for line in infile:
        el = line.strip().split()
        assert(len(el) == 3)
        hdf5_file_list[hg2na[el[2]]].append((el[0],el[1]))

#assert(len(hdf5_file_list[]) == num_hdf5_archives)

# for a single individual, merge all of their sorted, mapped pacbio bams into one big one.
rule merge_giab_pacbio_blasr:
    params: job_name = 'merge_giab_pacbio_blasr.{individual}',
    input: lambda wildcards: expand('data/{individual}/raw_pacbio/sorted_bam/archive{archive_number}.bam',individual=wildcards.individual, archive_number=list(range(num_hdf5_archives[wildcards.individual])))
    output: 'study/data/{individual}/aligned_reads/pacbio/pacbio.blasr.all.giab_full_coverage.bam'
    shell: '{SAMTOOLS} merge -@ 8 -O bam {output} {input}'

# sort a single pacbio bam
rule sort_giab_pacbio_blasr:
    params: job_name = 'sort_giab_pacbio_blasr.{individual}.archive{archive_number}',
    input: 'data/{individual}/raw_pacbio/aligned_bam/archive{archive_number}.bam'
    output: temp('data/{individual}/raw_pacbio/sorted_bam/archive{archive_number}.bam'),
    shell: '{SAMTOOLS} sort -T {output}.TMP -@ 4 -m 3G {input} > {output}'

# align the converted subread bam file to the reference genome
rule align_giab_pacbio_blasr:
    params: job_name = 'align_giab_pacbio_blasr.{individual}.archive{archive_number}'
    input:
        bam   = 'data/{individual}/raw_pacbio/unaligned_bam/archive{archive_number}.subreads.bam',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_sa = 'data/genomes/hs37d5.fa.sawriter.sa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
    output: bam = temp('data/{individual}/raw_pacbio/aligned_bam/archive{archive_number}.bam')
    shell: '{BLASR} {input.bam} {input.hs37d5} --sa {input.hs37d5_sa} --nproc 4 --bam --out {output}'

# convert the old hdf5 file format for pacbio reads to the new one compatible with
rule bax2bam:
    params: job_name = 'bax2bam.{individual}.archive{archive_number}',
            archive_dir = 'data/{individual}/raw_pacbio/hdf5/archive{archive_number}'
    input: 'data/{individual}/raw_pacbio/hdf5/archive{archive_number}/archive.tar.gz'
    output: temp('data/{individual}/raw_pacbio/unaligned_bam/archive{archive_number}.subreads.bam'),
            temp('data/{individual}/raw_pacbio/unaligned_bam/archive{archive_number}.subreads.bam.pbi'),
            temp('data/{individual}/raw_pacbio/unaligned_bam/archive{archive_number}.scraps.bam'),
            temp('data/{individual}/raw_pacbio/unaligned_bam/archive{archive_number}.scraps.bam.pbi')
    shell:
        '''
        cd {params.archive_dir}
        tar -xzf archive.tar.gz
        {BAX2BAM} *.bax.h5 -o ../../unaligned_bam/archive{wildcards.archive_number}
        rm *
        '''

rule download_pacbio_hdf5:
    params: job_name = 'download_pacbio_{individual}_hdf5.archive{archive_number}',
            archive_dir = 'data/{individual}/raw_pacbio/hdf5/archive{archive_number}'
    output: archive = 'data/{individual}/raw_pacbio/hdf5/archive{archive_number}/archive.tar.gz',
    run:
        a = int(wildcards.archive_number)
        hdf5_url = hdf5_file_list[wildcards.individual][a][0]
        hdf5_md5 = hdf5_file_list[wildcards.individual][a][1]
        shell('''
        wget {hdf5_url} -O {output.archive}
        md5sum -c <<<"{hdf5_md5} *{output.archive}"
        ''')
