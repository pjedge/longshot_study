import plot_vcfeval_precision_recall as plot_vcfeval

include: "simulation.snakefile"
include: "NA12878.snakefile"

# DATA URLs
PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
#PACBIO_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam.bai'
GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
HG19_URL     = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit'
HS37D5_URL     = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
HG19_SDF_URL   = 'https://s3.amazonaws.com/rtg-datasets/references/hg19.sdf.zip'
Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'
#Illumina_30x_BAI_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam.bai'

# PATHS TO TOOLS
TWOBITTOFASTA = 'twoBitToFa' # can be downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'
SAMTOOLS       = '/opt/biotools/samtools/1.3/bin/samtools' # v1.3
FASTQ_DUMP     = 'fastq-dump' # v2.5.2
REAPER         = '../target/release/reaper' # v0.1
RTGTOOLS       = '/home/pedge/installed/rtg-tools-3.8.4/rtg' # v3.8.4, https://www.realtimegenomics.com/products/rtg-tools
BGZIP = 'bgzip'
TABIX = 'tabix'
FREEBAYES      = '/home/pedge/git/freebayes/bin/freebayes'

# PARAMS
#chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX']  # ['chr20']  #

# DEFAULT

rule all:
    input:
        'data/plots/NA12878_prec_recall_chr20.png',
        #'data/plots/simulation_prec_recall_all.png'

# NOTE!!! we are filtering out indels but also MNPs which we may call as multiple SNVs
# therefore this isn't totally correct and it'd probably be better to use ROC with indels+SNVs VCF.
rule vcfeval_rtgtools:
    params: job_name = 'vcfeval_rtgtools.{dataset}.{calls_name}.{chrom}',
            region_arg = lambda wildcards: '--region={}'.format(wildcards.chrom) if wildcards.chrom != 'all' else ''
    input:  calls_vcf = 'data/{dataset}/variants/{calls_name}/{chrom}.vcf.gz',
            calls_ix = 'data/{dataset}/variants/{calls_name}/{chrom}.vcf.gz.tbi',
            ground_truth = 'data/{dataset}/variants/ground_truth/ground_truth.SNVs_ONLY.vcf.gz',
            ground_truth_ix = 'data/{dataset}/variants/ground_truth/ground_truth.SNVs_ONLY.vcf.gz.tbi',
            region_filter ='data/{dataset}/variants/ground_truth/region_filter.bed',
            hg19_sdf = 'data/genomes/hg19.sdf'
    output: done = 'data/{dataset}/vcfeval/{calls_name}/{chrom}.done'
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcfeval \
        {params.region_arg} \
        -c {input.calls_vcf} \
        -b {input.ground_truth} \
        -e {input.region_filter} \
        -t {input.hg19_sdf} \
        -o data/{wildcards.dataset}/vcfeval/{wildcards.calls_name}/{wildcards.chrom};
        cp data/{wildcards.dataset}/vcfeval/{wildcards.calls_name}/{wildcards.chrom}/done {output.done};
        '''

# NOTE!!! we are filtering out indels but also MNPs which we may call as multiple SNVs
# therefore this isn't totally correct and it'd probably be better to use ROC with indels+SNVs VCF.
rule rtg_filter_SNVs_ground_truth:
    params: job_name = 'rtg_filter_SNVs_ground_truth.{dataset}',
    input:  vcfgz = 'data/{dataset}/variants/ground_truth/ground_truth.vcf.gz'
    output: vcfgz = 'data/{dataset}/variants/ground_truth/ground_truth.SNVs_ONLY.vcf.gz',
            #tbi = 'data/{dataset}/variants/ground_truth/ground_truth.SNVs_ONLY.vcf.gz.tbi'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --snps-only -i {input.vcfgz} -o {output.vcfgz}'

from filter_SNVs import filter_SNVs
rule filter_illumina_SNVs:
    params: job_name = 'filter_SNVs_illumina.{dataset}.chr{chrnum}',
    input:  vcf = 'data/{dataset}/variants/illumina_{cov}x/chr{chrnum}.vcf'
    output: vcf = 'data/{dataset}/variants/illumina_{cov}x.filtered/chr{chrnum}.vcf'
    run:
        cov_filter = int(float(wildcards.cov)*1.75)
        filter_SNVs(input.vcf, output.vcf, cov_filter, density_count=10, density_len=500, density_qual=50)

rule combine_chrom:
    params: job_name = 'combine_chroms.{dataset}.{calls_name}',
    input: expand('data/{{dataset}}/variants/{{calls_name}}/{chrom}.vcf',chrom=chroms)
    output: 'data/{dataset}/variants/{calls_name}/all.vcf'
    shell:
        '''
        grep -P '^#' {input[0]} > {output}; # grep header
        cat {input} | grep -Pv '^#' >> {output}; # cat files, removing the headers.
        '''

####################################################################################################################
chunksize = int(1e6)

hg19_size_list = [('chr1', 249250621),
 ('chr2', 243199373),
 ('chr3', 198022430),
 ('chr4', 191154276),
 ('chr5', 180915260),
 ('chr6', 171115067),
 ('chr7', 159138663),
 ('chr8', 146364022),
 ('chr9', 141213431),
 ('chr10', 135534747),
 ('chr11', 135006516),
 ('chr12', 133851895),
 ('chr13', 115169878),
 ('chr14', 107349540),
 ('chr15', 102531392),
 ('chr16', 90354753),
 ('chr17', 81195210),
 ('chr18', 78077248),
 ('chr19', 59128983),
 ('chr20', 63025520),
 ('chr21', 48129895),
 ('chr22', 51304566),
 ('chrX', 155270560),
 ('chrY', 59373566)]

for chrom, chrlen in hg19_size_list:
    for start in range(1,chrlen+1,chunksize):
        end = start+chunksize-1 if start+chunksize-1 < chrlen else chrlen
        chunklist.append((chrom,start,end))

regions = ['{}.{}.{}'.format(chrom,start,stop) for chrom,start,stop in chunklist]
chr20_regions = [x for x in regions if x[:6] == 'chr20.']

rule combine_vcfs:
    params: job_name = 'combine_chroms.{dataset}.cov{cov}.chr{chrnum}'
    input: expand('data/{{dataset}}/variants/reaper_{{cov,\d+}}x.{{options}}/split_chrom/chr{r}.vcf',r=chr20_regions)
    output: 'data/{dataset}/variants/reaper_{cov,\d+}x.{options}/chr20.vcf',
    shell:
        '''
        grep -P "^#" {input[0]} > {output}
        cat {input} | grep -Pv "^#" > {output}
        '''

rule run_reaper:
    params: job_name = 'reaper.{dataset}.cov{cov}.chr{chrnum}.{start}.{stop}',
    input:  bam = 'data/{dataset}/aligned_reads/pacbio/pacbio.{cov}x.bam',
            bai = 'data/{dataset}/aligned_reads/pacbio/pacbio.{cov}x.bam.bai',
            ref    = 'data/genomes/hg19.fa',
            ref_ix = 'data/genomes/hg19.fa.fai'
    output: vcf = 'data/{dataset}/variants/reaper_{cov,\d+}x.{options}/split_chrom/chr{chrnum}.{start}.{stop}.vcf',
    run:
        options_str = wildcards.options.replace('_',' ')
        shell('{REAPER} -r chr{wildcards.chrnum}:{wildcards.start}-{wildcards.stop} {options_str} --bam {input.bam} --ref {input.ref} --out {output.vcf}')

####################################################################################################################

'''
rule run_reaper:
    params: job_name = 'reaper.{dataset}.cov{cov}.chr{chrnum}',
    input:  bam = 'data/{dataset}/aligned_reads/pacbio/pacbio.{cov}x.bam',
            bai = 'data/{dataset}/aligned_reads/pacbio/pacbio.{cov}x.bam.bai',
            ref    = 'data/genomes/hg19.fa',
            ref_ix = 'data/genomes/hg19.fa.fai'
    output: vcf = 'data/{dataset}/variants/reaper_{cov,\d+}x.{options}/chr{chrnum}.vcf',
    run:
        options_str = wildcards.options.replace('_',' ')
        shell('{REAPER} -r chr{wildcards.chrnum} {options_str} --bam {input.bam} --ref {input.ref} --out {output.vcf}')
'''
# Call 30x Illumina variants
rule call_variants_Illumina:
    params: job_name = 'call_illumina.{dataset}.{cov}x',
    input: bam = 'data/{dataset}/aligned_reads/illumina/illumina.{cov}x.bam',
            bai = 'data/{dataset}/aligned_reads/illumina/illumina.{cov}x.bam.bai',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5 = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai'
    output: vcf = 'data/{dataset}/variants/illumina_{cov}x/chr{chrnum}.vcf'
    run:
        if wildcards.dataset in ['NA12878']:
            shell('''
            {FREEBAYES} -f {input.hs37d5} \
            --standard-filters \
            --region {wildcards.chrnum} \
             --genotype-qualities \
             {input.bam} \
             | awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}' \
              > {output.vcf}
            ''')
        else:
            shell('''
            {FREEBAYES} -f {input.hg19} \
            --standard-filters \
            --region chr{wildcards.chrnum} \
             --genotype-qualities \
             {input.bam} \
              > {output.vcf}
            ''')

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

rule index_vcf:
    params: job_name = lambda wildcards: 'tabix_vcf.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.vcf.gz'
    output: '{x}.vcf.gz.tbi'
    shell:  '{TABIX} -p vcf {input}'

# bgzip vcf
rule bgzip_vcf_calls:
    params: job_name = 'bgzip_vcf_calls.{dataset}.{calls_name}.{chrom}'
    input:  'data/{dataset}/variants/{calls_name}/{chrom}.vcf'
    output: 'data/{dataset}/variants/{calls_name}/{chrom,(all|chr.+)}.vcf.gz'
    shell:  '{BGZIP} -c {input} > {output}'

# bgzip vcf
rule bgzip_ground_truth:
    params: job_name = 'bgzip_ground_truth.{dataset}'
    input:  'data/{dataset}/variants/ground_truth/ground_truth.vcf'
    output: 'data/{dataset}/variants/ground_truth/ground_truth.vcf.gz'
    shell:  '{BGZIP} -c {input} > {output}'

# gunzip fastq
rule gunzip_fastq:
    params: job_name = lambda wildcards: 'gunzip_fastq.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.fastq.gz'
    output: '{x}.fastq'
    shell:  'gunzip {input}'

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

# BWA index
rule bwa_index_fasta:
    params: job_name = lambda wildcards: 'bwa_index_fasta.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.fa'
    output: '{x}.fa.bwt'
    shell: '{BWA} index {input}'
