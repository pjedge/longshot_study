
NA24149_HG38_ILLUMINA_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.60x.1.bam'
NA24149_HG38_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz'
NA24149_HG38_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh38/HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'

rule plot_pr_curve_NA24149_hg38:
    params: job_name = 'plot_pr_curve_NA24149.hg38.{aligner}.{chrom}',
            title = None
    input: il30 = 'data/NA24149.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
           pb32 = 'data/NA24149.hg38/vcfeval/longshot.pacbio.{aligner}.32x._/{chrom}',
           il30_cov = 'data/NA24149.hg38/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           pb32_cov = 'data/NA24149.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.32x.bam.median_coverage'
    output: png = 'data/plots/NA24149.hg38.{aligner}.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.il30, input.pb32],
                         ['Freebayes, Illumina {}x'.format(parse_int_file(input.il30_cov)),
                         'Longshot, PacBio {}x'.format(parse_int_file(input.pb32_cov))],
                          output.png,params.title,
                          xlim=(0.6,1.0),
                          ylim=(0.975,1.0))

# DOWNLOAD 60x Illumina reads
rule download_Illumina_reads_NA24149_hg38:
    params: job_name = 'DOWNLOAD_Illumina_60x_NA24149.hg38',
    output: bam = 'data/NA24149.hg38/aligned_reads/illumina/illumina.aligned.all.60x.bam',
    shell: 'wget {NA24149_HG38_ILLUMINA_60x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24149_hg38:
    params: job_name = 'DOWNLOAD_GIAB_BED_NA24149.hg38',
    output: 'data/NA24149.hg38/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24149_HG38_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24149_hg38:
    params: job_name = 'DOWNLOAD_GIAB_VCF_NA24149.hg38'
    output: 'data/NA24149.hg38/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24149_HG38_GIAB_VCF_URL} -O {output}'

# RENAME PACBIO BAM
rule rename_pacbio_NA24149_hg38:
    params: job_name = 'rename_pacbio_NA24149.hg38.{aligner}'
    input: bam = 'data/NA24149.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.giab_full_coverage.bam'
    output: bam = 'data/NA24149.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.32x.bam',
    shell: 'mv {input.bam} {output.bam}'
