import statistics

NA12878_PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
NA12878_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
NA12878_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
NA12878_Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'

rule plot_pr_curve_NA12878:
    params: job_name = 'plot_pr_curve_NA12878.1000g.{aligner}.{chrom}',
            title = None
    input:
        il30 = 'data/NA12878.1000g/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
        pb30 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.30x._/{chrom}',
        pb44 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.44x._/{chrom}',
        il30_cov = 'data/NA12878.1000g/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
        pb30_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage',
        pb44_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.44x.bam.median_coverage'
    output:
        png = 'data/plots/NA12878.1000g.{aligner}.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.il30, input.pb30, input.pb44],
                         ['Freebayes, Illumina {}x'.format(parse_int_file(input.il30_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb30_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb44_cov))],
                          output.png,params.title,
                          colors=['r','#8080ff','#3333ff'],
                          xlim=(0.8,1.0),
                          ylim=(0.985,1.0))

rule compute_N50_haplotype_reads:
    params: job_name = 'compute_N50_haplotype_reads'
    input:  'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.{group}.bam'
    output: 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.{group,(hap1|hap2|unassigned)}.tlens.txt'
    shell: '{SAMTOOLS} view {input} | cut -f 9 > {output}'

def n50(lst):
    lst.sort(reverse=True)
    total = sum(lst)
    half_total = total/2.0
    count = 0
    for l in lst:
        count += l
        if count > half_total:
            return l
    return 0

def parse_file_lst(fname):
    lst = []
    with open(fname,'r') as inf:
        for line in inf:
            lst.append(int(line.strip()))
    return lst


rule get_haplotype_bam_stats:
    params: job_name = 'get_haplotype_bam_stats'
    input: hap1 = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.hap1.tlens.txt',
           hap2 = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.hap2.tlens.txt',
           una = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.unassigned.tlens.txt',
    output: txt = 'data/output/NA12878_chr1_45x_haplotype_assigned_N50_analysis.txt'
    run:
        hap_lst = parse_file_lst(input.hap1) + parse_file_lst(input.hap2)
        una_lst = parse_file_lst(input.una)
        with open(output.txt,'w') as outf:
            print('''Haplotype-assigned reads for NA12878 chr1 (45x reads):
N50: {}
mean: {}
median: {}
unassigned reads for NA12878 chr1 (45x reads):
N50: {}
mean: {}
median: {}
            '''.format(n50(hap_lst), statistics.mean(hap_lst), statistics.median(hap_lst),
                       n50(una_lst), statistics.mean(una_lst), statistics.median(una_lst)), file=outf)


rule haplotype_separation_chr1:
    params: job_name = 'haplotype_separation_chr1'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
           bai = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.bai',
           cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.median_coverage',
           ref = 'data/genomes/hg19.fa'
    output: vcf = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.for_haplotype_separation/1.vcf',
            h1_bam = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.hap1.bam',
            h2_bam = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.hap2.bam',
            una_bam = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.unassigned.bam'
    run:
        median_cov = parse_int_file(input.cov)
        min_cov = int(median_cov - 5*sqrt(median_cov))
        max_cov = int(median_cov + 5*sqrt(median_cov))
        if min_cov < 0:
            min_cov = 0
        shell('{LONGSHOT} -r chr1 -F -c {min_cov} -C {max_cov} -s NA12878 --bam {input.bam} --ref {input.ref} --out {output.vcf} -p data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x')

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads_NA12878:
    params: job_name = 'download_Illumina_30x_NA12878.1000g',
    output: bam = 'data/NA12878.1000g/aligned_reads/illumina/illumina.aligned.all.30x.bam',
    shell: 'wget {NA12878_Illumina_30x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA12878:
    params: job_name = 'download_giab_bed_NA12878.1000g',
    output: 'data/NA12878.1000g/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA12878_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA12878:
    params: job_name = 'download_giab_VCF_NA12878.1000g',
    output: 'data/NA12878.1000g/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA12878_GIAB_VCF_URL} -O {output}'

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA12878:
    params: job_name = 'split_bam_pacbio_NA12878.1000g.{chrom}'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.split_chroms/{chrom}.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} chr{wildcards.chrom} > {output.bam}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA12878:
    params: job_name = 'subsample_pacbio_NA12878.1000g'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.{cov,30}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 44.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA12878:
    params: job_name = 'download_pacbio_NA12878.1000g'
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    shell: 'wget {NA12878_PACBIO_BAM_URL} -O {output.bam}'
