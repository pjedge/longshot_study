#NA24385_PACBIO_NGMLR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.bam'
#NA24385_HG38_PACBIO_BLASR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg002_gr37_'
NA24385_HG38_Illumina_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam'
NA24385_HG38_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz'
NA24385_HG38_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
NA24385_HG38_PACBIO_BLASR_BAM_URL = 'http://www-rcf.usc.edu/~mchaisso/hg002.passthrough.bam'

rule plot_pr_curve_NA24385_hg38:
    params: job_name = 'plot_pr_curve_NA24385.hg38',
            title = 'Precision Recall Curve for Reaper on NA24385: PacBio Reads vs Standard Illumina'
    input:
        reaper20_rtg = 'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.20x.-z/{chrom}.done',
        reaper30_rtg = 'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.30x.-z/{chrom}.done',
        reaper40_rtg = 'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.40x.-z/{chrom}.done',
        reaper50_rtg = 'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.50x.-z/{chrom}.done',
        reaper69_rtg = 'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.69x.-z/{chrom}.done',
        illumina_rtg = 'data/NA24385.hg38/vcfeval/illumina_30x.filtered/{chrom}.done'
    output:
        png = 'data/plots/NA24385.hg38_prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval(['data/NA24385.hg38/vcfeval/illumina_30x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.20x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.30x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.40x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.50x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385.hg38/vcfeval/reaper.filtered.pacbio.blasr.69x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 30x',
                                   'Reaper, PacBio 20x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 40x',
                                   'Reaper, PacBio 50x',
                                   'Reaper, PacBio 69x',],
                                   output.png,params.title,
                                   colors=['#ff0707','#9999ff','#8080ff','#6666ff','#3333ff','b'],
                                   xlim=(0.6,1.0),ylim=(0.975,1.0))

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads_NA24385_hg38:
    params: job_name = 'download_Illumina_60x_NA24385.hg38',
    output: bam = 'data/NA24385.hg38/aligned_reads/illumina/illumina.60x.bam',
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
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom}.{cov,(20|30|40|50)}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 69.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA24385_hg38_BLASR:
    params: job_name = 'split_bam_pacbio_NA24385.hg38.{chrom}'
    input: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam',
           bai = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam.bai'
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom,(\d+)}.69x.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} chr{wildcards.chrom} > {output.bam}'

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24385_hg38:
    params: job_name = 'download_pacbio_NA24385.hg38'
    output: bam = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam',
    shell: 'wget {NA24385_HG38_PACBIO_BLASR_BAM_URL} -O {output.bam}'
