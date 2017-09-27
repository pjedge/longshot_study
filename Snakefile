
# FILE PATHS
TWOBITTOFASTA = 'twoBitToFa' # can be downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'
PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
PACBIO_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam.bai'
GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
DATA        = '/oasis/tscc/data/pedge/reaper'
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
RTGTOOLS       = '~/git/rtg-tools-3.8.4/rtg' # https://www.realtimegenomics.com/products/rtg-tools
HG19_SDF_URL   = 'https://s3.amazonaws.com/rtg-datasets/references/hg19.sdf.zip'
BGZIP = 'bgzip'
TABIX = 'tabix'

Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'
Illumina_30x_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai'
# PARAMS
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX']

# DEFAULT
rule all:

    input: 'data/accuracy_reports/reaper.chr1/done',
           'data/accuracy_reports/illumina_30x.chr1/done'

rule accuracy_stats_rtgtools:
    params: job_name = 'accuracy_stats.{calls_name}.{chrom}',
    input:  calls_vcf = 'data/variants/{calls_name}/{chrom}.vcf.gz',
            calls_ix = 'data/variants/{calls_name}/{chrom}.vcf.gz.tbi',
            giab_ref = 'data/variants/GIAB/GIAB.NA12878.vcf.gz',
            giab_ix = 'data/variants/GIAB/GIAB.NA12878.vcf.gz.tbi',
            bed_filter ='data/variants/GIAB/high_confidence.bed',
            hg19_sdf = 'data/genomes/hg19.sdf'
    output: report = 'data/accuracy_reports/{calls_name}.{chrom}/done'
    run:
        '''
        {RTGTOOLS} vcfeval \
        --region={chrom} \
        -c {input.calls_vcf} \
        -b {input.giab_ref} \
        -e {input.bed_filter} \
        -t {input.hg19_sdf} \
        --output=data/accuracy_reports/{wildcards.calls_name}.{wildcards.chrom} \
        '''

rule run_reaper:
    params: job_name = 'reaper.{chrom}',
    input:  bam_file = 'data/aligned_reads/pacbio/{chrom}.bam',
            bai_file = 'data/aligned_reads/pacbio/{chrom}.bai',
            ref      = HG19,
            ref_ix   = HG19 + '.fai'
    output: vcf = 'data/variants/reaper/{chrom}.vcf',
    shell: '{REAPER} -z --bam {input.bam} --ref {input.ref} --out {output.vcf} --max_cov 90'

#rule intersect_beds_giab_exons:
#    params: job_name = 'intersect_beds_giab_exons',
#    input: giab_conf = 'data/variants/GIAB/high_confidence.chr_prefix.bed',
#           exons     = 'data/genome_features/hg19_exons.bed'
#    output: bed_filter = 'data/genome_features/exons_intersect_giab_conf.bed'
#    shell: 'bedtools intersect -a {input.giab_conf} -b {input.exons} | sort -k 1,1 -k2,2n > {output.bed_filter}'

# Call 30x Illumina variants
rule call_variants_Illumina:
    params: job_name = 'Call_Illumina_30x',
    input: bam = 'data/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam',
            bai = 'data/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam.bai',
            fa = HS37D5,
            fai = HS37D5 + '.fai'
    output: vcf = 'data/variants/illumina_30x/{chr}.vcf'
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
    output: bam = 'data/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam',
            bai = 'data/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam.bai',
    shell:
        '''
        wget {Illumina_30x_BAM_URL} -O {output.bam};
        wget {Illumina_30x_BAI_URL} -O {output.bai};
        '''

# download hg19 reference, for the aligned pacbio reads
rule download_hg19:
    params: job_name = 'download_hg19',
            url      = HG19_URL
    output: 'data/genomes/hg19.fa'
    shell:
        '''
        wget {params.url} -O {output}.2bit
        {TWOBITTOFASTA} {output}.2bit {output}
        '''

# download hg19 reference, for the aligned pacbio reads
rule download_hg19_sdf:
    params: job_name = 'download_hg19_sdf',
    output: 'data/genomes/hg19.sdf'
    shell:
        '''
        wget {HG19_SDF_URL} -O {output}.zip;
        unzip {output}.zip -d data/genomes
        '''

rule download_HS37D5:
    params: job_name = 'download_hs37d',
    output: 'data/genomes/hs37d5.fa'
    shell: 'wget {HS37D5_URL} -O {output}.gz; gunzip -c {output}.gz'

rule extract_chrom_GIAB_VCF:
    params: job_name = 'extract_chr_GIAB_VCF.{chrnum}',
    input:  'data/variants/GIAB/GIAB.NA12878.vcf'
    output: 'data/variants/GIAB/chr{chrnum}.vcf'
    shell: '''cat {input} | grep -P '^{chrnum}\t' > {output}'''

rule add_chr_GIAB_VCF:
    params: job_name = 'add_chr_GIAB_VCF',
    input:  gz = 'data/variants/GIAB/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
    output: vcf = 'data/variants/GIAB/GIAB.NA12878.vcf'
    shell:
        '''
        gunzip -c {input.gz} | \
        awk '{{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }}' \
        > {output.vcf}
        '''

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed:
    params: job_name = 'DOWNLOAD_GIAB',
    output: 'data/variants/GIAB/high_confidence.bed'
    shell: '''wget -q -O- {GIAB_HIGH_CONF_URL} | awk '{{print "chr" $0;}}' > {output}'''

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF:
    params: job_name = 'DOWNLOAD_GIAB',
    output: 'data/variants/GIAB/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
    shell: 'wget {GIAB_VCF_URL} -O {output}'

# DOWNLOAD PACBIO BAM
rule download_pacbio:
    params: job_name = 'download_pacbio'
    output: bam = 'data/aligned_reads/pacbio/pacbio.bam',
            bai = 'data/aligned_reads/pacbio/pacbio.bam.bai'
    shell:
        '''
        wget {PACBIO_BAM_URL} -O {output.bam}
        wget {PACBIO_BAI_URL} -O {output.bai}
        '''

# bgzip vcf
rule index_vcf:
    params: job_name = lambda wildcards: 'tabix_vcf.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.vcf.gz'
    output: '{x}.vcf.gz.tbi'
    shell:  '{TABIX} -p vcf {input}'

# bgzip vcf
rule bgzip_vcf:
    params: job_name = lambda wildcards: 'bgzip_vcf.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.vcf'
    output: '{x}.vcf.gz'
    shell:  '{BGZIP} -c {input} > {output}'

# index fasta reference
rule index_fasta:
    params: job_name = lambda wildcards: 'index_fa.{}'.format(str(wildcards.x).replace("/", "."))
    input:  fa  = '{x}.fa'
    output: fai = '{x}.fa.fai'
    shell:  '{SAMTOOLS} faidx {input.fa}'

#rule index_bam:
#    params: job_name = lambda wildcards: 'index_bam.{}'.format(str(wildcards.x).replace("/", "."))
#    input:  bam = '{x}.bam'
#    output: bai = '{x}.bam.bai'
#    shell:  '{SAMTOOLS} index {input.bam} {output.bai}'
