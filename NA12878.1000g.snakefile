import statistics

NA12878_PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
NA12878_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
NA12878_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
NA12878_Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'
NA12878_PACBIO_NGMLR_BAM_URL = 'http://www.bio8.cs.hku.hk/clairvoyante/bamUsed/PacBio-HG001-hg19/na12878_pacbio_mts_ngmlr-0.2.3_mapped.rg.bam'

rule plot_pr_curve_NA12878_3method:
    params: job_name = 'plot_pr_curve_NA12878.1000g.ngmlr.all',
            title = None
    input:
        wh44 = 'data/NA12878.1000g/vcfeval/whatshap.pacbio.ngmlr.44x.dp_filtered/all',
        cl44 = 'data/NA12878.1000g/vcfeval/clairvoyante.pacbio.ngmlr.44x.unfiltered/all',
        ls44 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.ngmlr.44x._/all',
    output:
        png = 'data/plots/3method_NA12878.1000g.ngmlr.prec_recall_all.png'
    run:
        ptf.plot_vcfeval([input.wh44, input.cl44, input.ls44],
                         ['WhatsHap',
                          'Clairvoyante',
                          'Longshot'],
                          output.png,params.title,
                          colors=['orange','k','blue'],
                          xlim=(0.9,1.0),
                          ylim=(0.9,1.0))

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
    input: hap1 = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.hap1.bam.read_lengths',
           hap2 = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.hap2.bam.read_lengths',
           una = 'data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x.unassigned.bam.read_lengths',
    output: txt = 'data/output/NA12878_chr1_45x_haplotype_assigned_N50_analysis.txt'
    run:
        hap_lst = parse_file_lst(input.hap1) + parse_file_lst(input.hap2)
        una_lst = parse_file_lst(input.una)
        total_hap = sum(hap_lst)
        total_una = sum(una_lst)
        fraction_bases_hap_assigned = total_hap / (total_hap + total_una)
        with open(output.txt,'w') as outf:
            print('''Haplotype-assigned reads for NA12878 chr1 (45x reads):
N50: {}
mean: {}
median: {}
{} of total bases are assigned to a haplotype.
unassigned reads for NA12878 chr1 (45x reads):
N50: {}
mean: {}
median: {}
            '''.format(n50(hap_lst), statistics.mean(hap_lst), statistics.median(hap_lst), fraction_bases_hap_assigned,
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
        max_cov = int(median_cov + 5*sqrt(median_cov))
        shell('{LONGSHOT} -r chr1 -F -C {max_cov} -s NA12878 --bam {input.bam} --ref {input.ref} --out {output.vcf} -p data/NA12878.1000g/aligned_reads/pacbio/haplotype_separated.pacbio.blasr.chr1.44x')

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

# SPLIT BAM
rule split_bam_NA12878:
    params: job_name = 'split_bam_{tech}.{aligner}.{cov}x_NA12878.1000g.{chrom}'
    input: bam = 'data/NA12878.1000g/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam',
           bai = 'data/NA12878.1000g/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam.bai',
    output: bam = 'data/NA12878.1000g/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.split_chroms/{chrom}.bam',
    run:
        w_chrom = 'chr'+wildcards.chrom if wildcards.tech == 'pacbio' else wildcards.chrom
        shell('{SAMTOOLS} view -hb {input.bam} {w_chrom} > {output.bam}')

# FILTER PACBIO BAM
rule filter_pacbio_NA12878:
    params: job_name = 'filter_pacbio_NA12878.1000g'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr_filtered.all.44x.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} -F 3844 -q 30 chr20 > {output.bam}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA12878:
    params: job_name = 'subsample_pacbio_NA12878.1000g'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.{cov,30}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 44.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# DOWNLOAD PACBIO BAM
rule download_pacbio_ngmlr_NA12878:
    params: job_name = 'download_pacbio_ngmlr_NA12878.1000g'
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.ngmlr.all.44x.bam',
    shell: 'wget {NA12878_PACBIO_NGMLR_BAM_URL} -O {output.bam}'

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA12878:
    params: job_name = 'download_pacbio_NA12878.1000g'
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    shell: 'wget {NA12878_PACBIO_BAM_URL} -O {output.bam}'
