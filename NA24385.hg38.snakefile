
NA24385_HG38_Illumina_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam'
NA24385_HG38_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz'
NA24385_HG38_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
#NA24385_HG38_PACBIO_BLASR_BAM_URL = 'http://www-rcf.usc.edu/~mchaisso/hg002.passthrough.bam'

rule plot_pr_curve_NA24385_hg38:
    params: job_name = 'plot_pr_curve_NA24385.hg38.{aligner}.{chrom}',
            title = None
    input: il30 = 'data/NA24385.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
           pb30 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.30x._/{chrom}',
           pb40 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.40x._/{chrom}',
           pb50 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.50x._/{chrom}',
           pb69 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.69x._/{chrom}',
           il30_cov = 'data/NA24385.hg38/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           pb30_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage',
           pb40_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.40x.bam.median_coverage',
           pb50_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.50x.bam.median_coverage',
           pb69_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.69x.bam.median_coverage',
    output: png = 'data/plots/NA24385.hg38.{aligner}.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.il30, input.pb30, input.pb40, input.pb50, input.pb69],
                         ['Freebayes, Illumina {}x'.format(parse_int_file(input.il30_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb30_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb40_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb50_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb69_cov))],
                          output.png,
                          params.title,
                          colors=['#ff0707','#8080ff','#6666ff','#3333ff','b'],
                          xlim=(0.6,1.0),
                          ylim=(0.975,1.0))

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads_NA24385_hg38:
    params: job_name = 'download_Illumina_60x_NA24385.hg38',
    output: bam = 'data/NA24385.hg38/aligned_reads/illumina/illumina.aligned.all.60x.bam',
    shell: 'wget {NA24385_HG38_Illumina_60x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24385_hg38:
    params: job_name = 'download_giab_bed.hg38',
    output: 'data/NA24385.hg38/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24385_HG38_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24385_hg38:
    params: job_name = 'download_giab_vcf_NA24385.hg38'
    output: 'data/NA24385.hg38/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24385_HG38_GIAB_VCF_URL} -O {output}'

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA24385_hg38_BLASR:
    params: job_name = 'split_bam_pacbio_NA24385.hg38.{chrom}'
    input: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam',
           bai = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam.bai'
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.split_chroms/{chrom,(\d+)}.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} chr{wildcards.chrom} > {output.bam}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA24385_hg38:
    params: job_name = 'subsample_pacbio_NA24385.hg38.{aligner}.{chrom}'
    input: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.{chrom}.69x.bam',
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.{chrom,(\d+|all)}.{cov,(20|30|40|50)}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 69.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# DOWNLOAD PACBIO BAM
#rule download_pacbio_NA24385_hg38:
#    params: job_name = 'download_pacbio_NA24385.hg38'
#    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam',
#    shell: 'wget {NA24385_HG38_PACBIO_BLASR_BAM_URL} -O {output.bam}'

# RENAME PACBIO BAM
rule rename_pacbio_NA24385_hg38:
    params: job_name = 'rename_pacbio_NA24385.hg38.blasr'
    input: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.giab_full_coverage.bam'
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam',
    shell: 'mv {input.bam} {output.bam}'
