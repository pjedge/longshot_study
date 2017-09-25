from count_accuracy import count_accuracy

# FILE PATHS
TWOBITTOFASTA = 'twoBitToFa' # can be downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'
PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
PACBIO_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam.bai'
GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
DATA        = '/oasis/tscc/scratch/pedge/reaper'
BWA            = '/home/pedge/installed/bwa'
SAMTOOLS       = '/home/pedge/git/samtools/samtools'
PICARD         = '/home/pedge/installed/picard/broadinstitute-picard-2a49ee2/dist/picard.jar'
bamtools       = '~/installed/bamtools/bin/bamtools'
FASTQ_DUMP     = 'fastq-dump'
#REF_URL        = 'ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
#REF_IX_URL     = 'ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai'
REAPER         = '../target/release/reaper'
#HG19           = '/oasis/tscc/scratch/pedge/data/genomes/hg19/hg19.fa'
#HG19_index     = '/oasis/tscc/scratch/pedge/data/genomes/hg19/hg19.fa.bwt'
HG19            = 'data/genomes/hg19.fa'
HG19_URL     = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit'
HS37D5         = 'data/genomes/hs37d5.fa'
HS37D5_URL     = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
#PYPY           = '/home/pedge/installed/pypy3.3-5.5-alpha-20161013-linux_x86_64-portable/bin/pypy3.3'


#COVERAGE_CUT = 63
Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'
Illumina_30x_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai'
# PARAMS
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX']

# DEFAULT
rule all:
    input: 'scratch/accuracy_reports/chr20.txt'

rule accuracy_stats:
    params: job_name = 'accuracy_stats.{chrom}',
    input:  reaper_vcf = '{DATA}/variants/reaper/{chrom}.vcf',
            illumina_vcf = '{DATA}/variants/illumina_30x/{chrom}.vcf',
            giab_ref = '{DATA}/variants/GIAB/{chrom}.vcf',
            bed_filter ='{DATA}/variants/GIAB/high_confidence.bed',
            #'{DATA}/variants/GIAB/high_confidence.chr_prefix.bed',
            hg19_fasta = HG19
    output: report = '{DATA}/accuracy_reports/{chrom}.txt'
    run:
        count_accuracy(input.reaper_vcf, input.illumina_vcf, input.giab_ref, input.bed_filter, wildcards.chrom, input.hg19_fasta, output.report)

rule run_reaper:
    params: job_name = 'reaper.{chrom}',
    input:  bam_file = '{DATA}/aligned_reads/pacbio/{chrom}.bam',
            bai_file = '{DATA}/aligned_reads/pacbio/{chrom}.bai',
            ref      = HG19
            ref_ix   = HG19 + '.fai'
    output: vcf = '{DATA}/variants/reaper/{chrom}.vcf',
    shell: '{REAPER} -z --bam {input.bam} --ref {input.ref} --out {output.vcf}'

#rule intersect_beds_giab_exons:
#    params: job_name = 'intersect_beds_giab_exons',
#    input: giab_conf = '{DATA}/variants/GIAB/high_confidence.chr_prefix.bed',
#           exons     = '{DATA}/genome_features/hg19_exons.bed'
#    output: bed_filter = '{DATA}/genome_features/exons_intersect_giab_conf.bed'
#    shell: 'bedtools intersect -a {input.giab_conf} -b {input.exons} | sort -k 1,1 -k2,2n > {output.bed_filter}'

# Call 30x Illumina variants
rule call_variants_Illumina:
    params: job_name = 'Call_Illumina_30x',
    input: bam = '{DATA}/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam',
            bai = '{DATA}/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam.bai',
            fa = HS37D5
            fai = HS37D5 + '.fai'
    output: vcf = '{DATA}/variants/illumina_30x/{chr}.vcf'
    shell:
        '''
        freebayes -f {HS37D5} \
        --standard-filters \
        --region {wildcards.chr} \
         {input.bam} \
         | grep -Pv '^#' \
         | sed -e "s/^/chr/" > {output.vcf}
        '''

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads:
    params: job_name = 'DOWNLOAD_Illumina_30x',
    output: bam = '{DATA}/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam',
            bai = '{DATA}/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam.bai',
    shell:
        '''
        wget {Illumina_30x_BAM_URL} -O {output.bam};
        wget {Illumina_30x_BAI_URL} -O {output.bai};
        '''

# download hg19 reference, for the aligned pacbio reads
rule download_hg19_pacbio:
    params: job_name = 'download_hg19',
            url      = HG19_URL
    output: 'data/genomes/hg19.fa'
    shell:
        '''
        wget {params.url} -O {output}.2bit
        {TWOBITTOFASTA} {output}.2bit {output}
        '''

rule download_HS37D5:
    output: '{DATA}/genomes/hs37d5.fa'
    shell: 'wget {HS37D5_URL} -O {output}.gz; gunzip -c {output}.gz'

# ADD CHR TO TG VARIANTS
rule extract_chrom_GIAB_VCF:
    params: job_name = 'extract_chr_GIAB_VCF.{chrnum}',
    input:  gz = '{DATA}/variants/GIAB/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
    output: vcf = '{DATA}/variants/GIAB/chr{chrnum}.vcf'
    shell: '''gunzip -c {input.gz} | grep -P '^{chrnum}\t' | awk '{print "chr" $0;}' > {output.vcf}'''

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed:
    params: job_name = 'DOWNLOAD_GIAB',
    output: '{DATA}/variants/GIAB/high_confidence.bed'
    shell: 'wget -q -O- {GIAB_HIGH_CONF_URL} | awk '{print "chr" $0;}' > {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF:
    params: job_name = 'DOWNLOAD_GIAB',
    output: '{DATA}/variants/GIAB/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
    shell: 'wget {GIAB_VCF_URL} -O {output}'

# DOWNLOAD PACBIO BAM
rule download_pacbio:
    params: job_name = 'download_pacbio'
    output: bam = '{DATA}/aligned_reads/pacbio/pacbio.bam',
            bai = '{DATA}/aligned_reads/pacbio/pacbio.bam.bai'
    shell:
        '''
        wget {PACBIO_BAM_URL} -O {output.bam}
        wget {PACBIO_BAI_URL} -O {output.bai}
        '''

# index fasta reference
rule index_fasta:
    params: job_name = lambda wildcards: 'index_fa.{}'.format(str(wildcards.x).replace("/", "."))
    input:  fa  = '{x}.fa'
    output: fai = '{x}.fa.fai'
    shell:  '{SAMTOOLS} faidx {input.fa}'

rule index_bam:
    params: job_name = lambda wildcards: 'index_bam.{}'.format(str(wildcards.x).replace("/", "."))
    input:  bam = '{x}.bam'
    output: bai = '{x}.bam.bai'
    shell:  '{SAMTOOLS} index {input.bam} {output.bai}'
