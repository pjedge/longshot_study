
NA24143_HG38_ILLUMINA_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.60x.1.bam'
NA24143_HG38_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf.vcf.gz'
NA24143_HG38_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh38/HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
NA24143_HG38_PACBIO_BLASR_BAM_URL = None #'http://www-rcf.usc.edu/~mchaisso/hg002.passthrough.bam'

rule plot_pr_curve_NA24143_hg38:
    params: job_name = 'plot_pr_curve_NA24143.hg38',
            title = 'Precision Recall Curve for Reaper on NA24143: PacBio Reads vs Standard Illumina'
    input:
        reaper30_rtg = 'data/NA24143.hg38/vcfeval/reaper.pacbio.blasr.30x.-z/{chrom}',
        illumina_rtg = 'data/NA24143.hg38/vcfeval/illumina_30x.filtered/{chrom}'
    output:
        png = 'data/plots/NA24143.hg38_prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval(['data/NA24143.hg38/vcfeval/illumina_30x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA24143.hg38/vcfeval/reaper.pacbio.blasr.30x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 30x',
                                   'Reaper, PacBio 30x'],
                                   output.png,params.title,
                                   xlim=(0.6,1.0),ylim=(0.975,1.0))

# DOWNLOAD 60x Illumina reads
rule download_Illumina_reads_NA24143_hg38:
    params: job_name = 'DOWNLOAD_Illumina_60x_NA24143.hg38',
    output: bam = 'data/NA24143.hg38/aligned_reads/illumina/illumina.60x.bam',
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

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA24143_blasr_hg38:
    params: job_name = 'split_bam_pacbio_NA24143.hg38.{chrom}'
    input: bam = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.all.30x.bam',
           bai = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.bai'
    output: bam = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom,(\d+)}.30x.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} chr{wildcards.chrom} > {output.bam}'

# DOWNLOAD PACBIO BAM
#rule download_pacbio_NA24143_hg38:
#    params: job_name = 'download_pacbio_NA24143.hg38'
#    output: bam = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.all.30x.bam',
#    shell: 'wget {NA24143_HG38_PACBIO_BLASR_BAM_URL} -O {output.bam}'

# RENAME PACBIO BAM
rule rename_pacbio_NA24143_hg38:
    params: job_name = 'rename_pacbio_NA24143.hg38'
    input: bam = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.all.giab_full_coverage.bam'
    output: bam = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.all.30x.bam',
    shell: 'mv {input.bam} {output.bam}'
