
NA24385_HG38_Illumina_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam'
NA24385_HG38_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz'
NA24385_HG38_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
NA24385_HG38_PACBIO_BLASR_BAM_URL = 'http://www-rcf.usc.edu/~mchaisso/hg002.passthrough.bam'

rule plot_pr_curve_NA24385_hg38:
    params: job_name = 'plot_pr_curve_NA24385.hg38.{aligner}.{chrom}',
            title = 'Precision Recall Curve for Reaper on NA24385: PacBio Reads vs Standard Illumina'
    input: il30 = 'data/NA24385.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}'
           pb20 = 'data/NA24385.hg38/vcfeval/reaper.pacbio.blasr.20x._/{chrom}',
           pb30 = 'data/NA24385.hg38/vcfeval/reaper.pacbio.blasr.30x._/{chrom}',
           pb40 = 'data/NA24385.hg38/vcfeval/reaper.pacbio.blasr.40x._/{chrom}',
           pb50 = 'data/NA24385.hg38/vcfeval/reaper.pacbio.blasr.50x._/{chrom}',
           pb69 = 'data/NA24385.hg38/vcfeval/reaper.pacbio.blasr.69x._/{chrom}',
           il30_cov = 'data/NA24385.hg38/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           pb20_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.20x.bam.median_coverage',
           pb30_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage',
           pb40_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.40x.bam.median_coverage',
           pb50_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.50x.bam.median_coverage',
           pb69_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.69x.bam.median_coverage',
    output: png = 'data/plots/NA24385.hg38.{aligner}.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.il30, input.pb20, input.pb30, input.pb40, input.pb50, input.pb69],
                         ['Freebayes, Illumina {}x'.format(parse_int_file(input.il30_cov)),
                          'Reaper, PacBio {}x'.format(parse_int_file(input.pb20_cov)),
                          'Reaper, PacBio {}x'.format(parse_int_file(input.pb30_cov)),
                          'Reaper, PacBio {}x'.format(parse_int_file(input.pb40_cov)),
                          'Reaper, PacBio {}x'.format(parse_int_file(input.pb50_cov)),
                          'Reaper, PacBio {}x'.format(parse_int_file(input.pb69_cov))],
                          output.png,
                          params.title,
                          colors=['#ff0707','#9999ff','#8080ff','#6666ff','#3333ff','b'],
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

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA24385_hg38:
    params: job_name = 'subsample_pacbio_NA24385.hg38.{chrom}'
    input: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom}.69x.bam',
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom,\d+}.{cov,(20|30|40|50)}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 69.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24385_hg38:
    params: job_name = 'download_pacbio_NA24385.hg38'
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam',
    shell: 'wget {NA24385_HG38_PACBIO_BLASR_BAM_URL} -O {output.bam}'
