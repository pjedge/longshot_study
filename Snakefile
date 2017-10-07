
# DATA URLs
PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
PACBIO_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam.bai'
GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
HG19_URL     = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit'
HS37D5_URL     = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
HG19_SDF_URL   = 'https://s3.amazonaws.com/rtg-datasets/references/hg19.sdf.zip'
Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'
Illumina_30x_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai'

# PATHS TO TOOLS
TWOBITTOFASTA = 'twoBitToFa' # can be downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'
SAMTOOLS       = '/opt/biotools/samtools/1.3/bin/samtools' # v1.3
FASTQ_DUMP     = 'fastq-dump' # v2.5.2
REAPER         = '../target/release/reaper' # v0.1
RTGTOOLS       = '/home/pedge/installed/rtg-tools-3.8.4/rtg' # v3.8.4, https://www.realtimegenomics.com/products/rtg-tools
BGZIP = 'bgzip'
TABIX = 'tabix'

# PARAMS
chroms = ['chr20']  #['chr{}'.format(i) for i in range(1,23)] + ['chrX']

# DEFAULT
methods = [
'reaper_-x_-n_-m_1_-W_500_-B_30_-i_-C_77', # max alignment, no short haps, no long haps
'reaper_-x_-n_-W_500_-B_30_-i_-C_77', # max alignment, no long haps
'reaper_-z_-n_-W_500_-B_30_-i_-C_77', # all alignment, no long haps
'reaper_-z_-W_500_-B_30_-i_-C_77', # all alignment
'illumina_30x.filtered'
]

rule all:
    input: expand('data/vcfeval/{m}.{c}.done',m=methods,c=['chr20'])
    #expand('data/vcfeval/{m}.{c}.done',m=methods,c=chroms+['all'])

#rule vcfeval_rtgtools_all:
#    params: job_name = 'vcfeval_rtgtools.{calls_name}.all',
#    input:  calls_vcf = 'data/variants/{calls_name}/all.vcf.gz',
#            calls_ix = 'data/variants/{calls_name}/all.vcf.gz.tbi',
#            giab_ref = 'data/variants/GIAB/NA12878.vcf.gz',
#            giab_ix = 'data/variants/GIAB/NA12878.vcf.gz.tbi',
#            bed_filter ='data/variants/GIAB/high_confidence.bed',
#            hg19_sdf = 'data/genomes/hg19.sdf'
#    output: done = 'data/vcfeval/{calls_name}.all.done'
#    shell:
#        '''
#        {RTGTOOLS} RTG_MEM=12g vcfeval \
#        --region={wildcards.chrom} \
#        -c {input.calls_vcf} \
#        -b {input.giab_ref} \
#        -e {input.bed_filter} \
#        -t {input.hg19_sdf} \
#        -o data/vcfeval/{wildcards.calls_name}.{wildcards.chrom} \
#        --output-mode=roc-only;
#        cp data/vcfeval/{wildcards.calls_name}.{wildcards.chrom}/done {output.done};
#        '''


# NOTE!!! we are filtering out indels but also MNPs which we may call as multiple SNVs
# therefore this isn't totally correct and it'd probably be better to use ROC with indels+SNVs VCF.
rule vcfeval_rtgtools:
    params: job_name = 'vcfeval_rtgtools.{calls_name}.{chrom}',
            region_arg = lambda wildcards: '--region={}'.format(wildcards.chrom) if wildcards.chrom != 'all' else ''
    input:  calls_vcf = 'data/variants/{calls_name}/{chrom}.vcf.gz',
            calls_ix = 'data/variants/{calls_name}/{chrom}.vcf.gz.tbi',
            giab_ref = 'data/variants/GIAB/NA12878.SNVs_ONLY.vcf.gz',
            giab_ix = 'data/variants/GIAB/NA12878.SNVs_ONLY.vcf.gz.tbi',
            bed_filter ='data/variants/GIAB/high_confidence.bed',
            hg19_sdf = 'data/genomes/hg19.sdf'
    output: done = 'data/vcfeval/{calls_name}.{chrom}.done'
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcfeval \
        {params.region_arg} \
        -c {input.calls_vcf} \
        -b {input.giab_ref} \
        -e {input.bed_filter} \
        -t {input.hg19_sdf} \
        -o data/vcfeval/{wildcards.calls_name}.{wildcards.chrom};
        cp data/vcfeval/{wildcards.calls_name}.{wildcards.chrom}/done {output.done};
        '''

# NOTE!!! we are filtering out indels but also MNPs which we may call as multiple SNVs
# therefore this isn't totally correct and it'd probably be better to use ROC with indels+SNVs VCF.
rule rtg_filter_SNVs_GIAB:
    params: job_name = 'rtg_filter_SNVs_GIAB',
    input:  vcfgz = 'data/variants/GIAB/NA12878.vcf.gz'
    output: vcfgz = 'data/variants/GIAB/NA12878.SNVs_ONLY.vcf.gz',
            tbi = 'data/variants/GIAB/NA12878.SNVs_ONLY.vcf.gz.tbi'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --snps-only -i {input.vcfgz} -o {output.vcfgz}'

from filter_SNVs import filter_SNVs
rule filter_illumina_SNVs:
    params: job_name = 'filter_SNVs_illumina.{chrom}',
    input:  vcf = 'data/variants/illumina_30x/{chrom}.vcf'
    output: vcf = 'data/variants/illumina_30x.filtered/{chrom}.vcf'
    run:
        filter_SNVs(input.vcf, output.vcf, 52, density_count=10, density_len=500, density_qual=50)

rule combine_chrom:
    params: job_name = 'combine_chroms.{calls_name}',
    input: expand('data/variants/{{calls_name}}/{chrom}.vcf',chrom=chroms)
    output: 'data/variants/{calls_name}/all.vcf'
    shell:
        '''
        grep -P '^#' {input[0]} > {output}; # grep header
        cat {input} | grep -Pv '^#' >> {output}; # cat files, removing the headers.
        '''

rule run_reaper:
    params: job_name = 'reaper.chr{chrnum}',
    input:  bam = 'data/aligned_reads/pacbio/pacbio.bam',
            bai = 'data/aligned_reads/pacbio/pacbio.bam.bai',
            ref    = 'data/genomes/hg19.fa',
            ref_ix = 'data/genomes/hg19.fa.fai'
    output: vcf = 'data/variants/reaper_{options}/chr{chrnum}.vcf',
    run:
        options_str = wildcards.options.replace('_',' ')
        shell('{REAPER} -r chr{wildcards.chrnum} {options_str} --bam {input.bam} --ref {input.ref} --out {output.vcf}')

#rule intersect_beds_giab_exons:
#    params: job_name = 'intersect_beds_giab_exons',
#    input: giab_conf = 'data/variants/GIAB/high_confidence.chr_prefix.bed',
#           exons     = 'data/genome_features/hg19_exons.bed'
#    output: bed_filter = 'data/genome_features/exons_intersect_giab_conf.bed'
#    shell: 'bedtools intersect -a {input.giab_conf} -b {input.exons} | sort -k 1,1 -k2,2n > {output.bed_filter}'

# Call 30x Illumina variants
rule call_variants_Illumina:
    params: job_name = 'call_illumina_30x',
    input: bam = 'data/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam',
            bai = 'data/aligned_reads/illumina_30x/NIST_NA12878_HG001_HiSeq_300x_RMNISTHS_30xdownsample.bam.bai',
            fa = 'data/genomes/hs37d5.fa',
            fai = 'data/genomes/hs37d5.fa.fai'
    output: vcf = 'data/variants/illumina_30x/chr{chrnum}.vcf'
    shell:
        '''
        freebayes -f {input.fa} \
        --standard-filters \
        --region {wildcards.chrnum} \
         --genotype-qualities \
         {input.bam} \
         | awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}' \
          > {output.vcf}
        '''
        #--max-coverage 60 \

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
    shell: 'wget {HS37D5_URL} -O {output}.gz; gunzip {output}.gz'

rule extract_chrom_GIAB_VCF:
    params: job_name = 'extract_chr_GIAB_VCF.{chrnum}',
    input:  'data/variants/GIAB/NA12878.vcf'
    output: 'data/variants/GIAB/chr{chrnum}.vcf'
    shell: '''cat {input} | grep -P '^{chrnum}\t' > {output}'''

rule add_chr_GIAB_VCF:
    params: job_name = 'add_chr_GIAB_VCF',
    input:  gz = 'data/variants/GIAB/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
    output: vcf = 'data/variants/GIAB/NA12878.vcf'
    shell:
        '''
        gunzip -c {input.gz} | \
        awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}' \
        > {output.vcf}
        '''

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed:
    params: job_name = 'DOWNLOAD_GIAB_BED',
    output: 'data/variants/GIAB/high_confidence.bed'
    shell: '''wget -q -O- {GIAB_HIGH_CONF_URL} | awk '{{print "chr" $0;}}' > {output}'''

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF:
    params: job_name = 'DOWNLOAD_GIAB_VCF',
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
