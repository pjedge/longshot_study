#NA24385_PACBIO_NGMLR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.bam'
NA24385_PACBIO_BLASR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg002_gr37_'
NA24385_Illumina_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam'
NA24385_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz'
NA24385_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'

rule plot_pr_curve_NA24385:
    params: job_name = 'plot_pr_curve_NA24385',
            title = 'Precision Recall Curve for Reaper on NA24385: PacBio Reads vs Standard Illumina'
    input:
        reaper20_rtg = 'data/NA24385/vcfeval/reaper.pacbio.blasr.20x.-z/{chrom}.done',
        reaper30_rtg = 'data/NA24385/vcfeval/reaper.pacbio.blasr.30x.-z/{chrom}.done',
        reaper40_rtg = 'data/NA24385/vcfeval/reaper.pacbio.blasr.40x.-z/{chrom}.done',
        reaper50_rtg = 'data/NA24385/vcfeval/reaper.pacbio.blasr.50x.-z/{chrom}.done',
        reaper69_rtg = 'data/NA24385/vcfeval/reaper.pacbio.blasr.69x.-z/{chrom}.done',
        illumina_rtg = 'data/NA24385/vcfeval/illumina_60x.filtered/{chrom}.done'
    output:
        png = 'data/plots/NA24385_prec_recall_{chrom}.png'
    run:
        plot_vcfeval.plot_vcfeval(['data/NA24385/vcfeval/illumina_60x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA24385/vcfeval/reaper.pacbio.blasr.20x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385/vcfeval/reaper.pacbio.blasr.30x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385/vcfeval/reaper.pacbio.blasr.40x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385/vcfeval/reaper.pacbio.blasr.50x.-z/{}'.format(wildcards.chrom),
                                   'data/NA24385/vcfeval/reaper.pacbio.blasr.69x.-z/{}'.format(wildcards.chrom),],
                                   ['Freebayes, Illumina 60x',
                                   'Reaper, PacBio 20x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 40x',
                                   'Reaper, PacBio 50x',
                                   'Reaper, PacBio 69x',],
                                   output.png,params.title)

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads_NA24385:
    params: job_name = 'DOWNLOAD_Illumina_60x_NA24385',
    output: bam = 'data/NA24385/aligned_reads/illumina/illumina.60x.bam',
    shell: 'wget {NA24385_Illumina_60x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24385:
    params: job_name = 'DOWNLOAD_GIAB_BED_NA24385',
    output: 'data/NA24385/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24385_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24385:
    params: job_name = 'DOWNLOAD_GIAB_VCF_NA24385'
    output: 'data/NA24385/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24385_GIAB_VCF_URL} -O {output}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA24385:
    params: job_name = 'subsample_pacbio_NA24385.{chrom}'
    input: bam = 'data/NA24385/aligned_reads/pacbio/pacbio.blasr.{chrom}.69x.bam',
    output: bam = 'data/NA24385/aligned_reads/pacbio/pacbio.blasr.{chrom}.{cov}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 69.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24385:
    params: job_name = 'download_pacbio_NA24385.{chrom}'
    output: bam = 'data/NA24385/aligned_reads/pacbio/pacbio.blasr.{chrom}.69x.bam',
    shell: 'wget {NA24385_PACBIO_BLASR_BAM_URL}{wildcards.chrom}.bam -O {output.bam}'

# DOWNLOAD PACBIO BAM
#rule download_pacbio_NA24385:
#    params: job_name = 'download_pacbio_NA24385'
#    output: bam = 'data/NA24385/aligned_reads/pacbio/pacbio.69x.bam',
#    shell: 'wget {NA24385_PACBIO_BAM_URL} -O {output.bam}'
