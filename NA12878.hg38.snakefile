
NA12878_HG38_ONT_BAM_URL = 'https://s3.amazonaws.com/nanopore-human-wgs/rel5-guppy-0.3.0-chunk10k.sorted.bam'
NA12878_HG38_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed'
NA12878_HG38_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'

rule plot_pr_curve_NA12878_pacbio_vs_ONT:
    params: job_name = 'plot_pr_curve_NA12878_pacbio_vs_ONT',
            title = None
    input:
        pb30 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.30x._/all',
        pb30_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
        ont30 = 'data/NA12878.hg38/vcfeval/longshot.ont.minimap2.30x._/all',
        ont30_cov = 'data/NA12878.hg38/aligned_reads/ont/ont.minimap2.all.30x.bam.median_coverage',
    output:
        png = 'data/plots/pacbio_vs_ont_PR_curve_all.png'
    run:
        ptf.plot_vcfeval([input.pb30, input.ont30],
                         ['Longshot, PacBio {}x'.format(parse_int_file(input.pb30_cov)),
                          'Longshot, Oxford Nanopore {}x'.format(parse_int_file(input.pb30_cov))],
                          output.png,params.title,
                          colors=['b','g'],
                          xlim=(0,1.0),
                          ylim=(0,1.0))

rule plot_pr_curve_NA12878_ONT:
    params: job_name = 'plot_pr_curve_NA12878_ONT',
            title = None
    input:
        ont30 = 'data/NA12878.hg38/vcfeval/longshot.ont.minimap2.30x._/all',
        ont30_cov = 'data/NA12878.hg38/aligned_reads/ont/ont.minimap2.all.30x.bam.median_coverage',
    output:
        png = 'data/plots/NA12878.hg38_ONT_PR_curve_all.png'
    run:
        ptf.plot_vcfeval([input.ont30],
                         ['Longshot, Oxford Nanopore {}x'.format(parse_int_file(input.pb30_cov))],
                          output.png, params.title,
                          colors=['g'],
                          xlim=(0,1.0),
                          ylim=(0,1.0))

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA12878_hg38:
    params: job_name = 'download_giab_bed_NA12878.hg38',
    output: 'data/NA12878.hg38/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA12878_HG38_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA12878_hg38:
    params: job_name = 'download_giab_VCF_NA12878.hg38',
    output: 'data/NA12878.hg38/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA12878_HG38_GIAB_VCF_URL} -O {output}'

# SPLIT BAM
rule split_bam_NA12878_hg38:
    params: job_name = 'split_bam_{tech}.{aligner}.{cov}x_NA12878.hg38.{chrom}'
    input: bam = 'data/NA12878.hg38/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam',
    output: bam = 'data/NA12878.hg38/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.split_chroms/{chrom}.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} {w_chrom} > {output.bam}'

# DOWNLOAD ONT BAM
rule download_ont_NA12878_hg38:
    params: job_name = 'download_ont_NA12878.hg38'
    output: bam = 'data/NA12878.hg38/aligned_reads/ont/ont.minimap2.all.30x.bam',
    shell: 'wget {NA12878_HG38_ONT_BAM_URL} -O {output.bam}'
