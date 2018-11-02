
NA24143_HG38_ILLUMINA_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.60x.1.bam'
NA24143_HG38_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz'
NA24143_HG38_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'

rule plot_pr_curve_NA24143_hg38_Q30:
    params: job_name = 'plot_pr_curve_Q30_NA24143.hg38.{aligner}.{chrom}',
            title = None
    input: il30 = 'data/NA24143.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
           pb30 = 'data/NA24143.hg38/vcfeval/longshot.pacbio.{aligner}.30x.-Q_30/{chrom}',
           il30_cov = 'data/NA24143.hg38/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           pb30_cov = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage'
    output: png = 'data/plots/Q30_NA24143.hg38.{aligner}.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.il30, input.pb30],
                         ['Freebayes, Illumina {}x'.format(parse_int_file(input.il30_cov)),
                         'Longshot, PacBio {}x'.format(parse_int_file(input.pb30_cov))],
                          output.png,params.title,
                          xlim=(0.6,1.0),
                          ylim=(0.975,1.0))

rule plot_pr_curve_NA24143_hg38:
    params: job_name = 'plot_pr_curve_NA24143.hg38.{aligner}.{chrom}',
            title = None
    input: il30 = 'data/NA24143.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
           pb30 = 'data/NA24143.hg38/vcfeval/longshot.pacbio.{aligner}.30x._/{chrom}',
           il30_cov = 'data/NA24143.hg38/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           pb30_cov = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage'
    output: png = 'data/plots/NA24143.hg38.{aligner}.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.il30, input.pb30],
                         ['Freebayes, Illumina {}x'.format(parse_int_file(input.il30_cov)),
                         'Longshot, PacBio {}x'.format(parse_int_file(input.pb30_cov))],
                          output.png,params.title,
                          xlim=(0.6,1.0),
                          ylim=(0.975,1.0))

# DOWNLOAD 60x Illumina reads
rule download_Illumina_reads_NA24143_hg38:
    params: job_name = 'DOWNLOAD_Illumina_60x_NA24143.hg38',
    output: bam = 'data/NA24143.hg38/aligned_reads/illumina/illumina.aligned.all.60x.bam',
    shell: 'wget {NA24143_HG38_ILLUMINA_60x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24143_hg38:
    params: job_name = 'DOWNLOAD_GIAB_BED_NA24143.hg38',
    output: 'data/NA24143.hg38/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24143_HG38_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24143_hg38:
    params: job_name = 'DOWNLOAD_GIAB_VCF_NA24143.hg38'
    output: 'data/NA24143.hg38/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24143_HG38_GIAB_VCF_URL} -O {output}'

# RENAME PACBIO BAM
rule rename_pacbio_NA24143_hg38:
    params: job_name = 'rename_pacbio_NA24143.hg38.{aligner}'
    input: bam = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.giab_full_coverage.bam'
    output: bam = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam',
    shell: 'mv {input.bam} {output.bam}'
