
NA24149_1000G_ILLUMINA_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.60x.1.bam'
NA24149_1000G_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'
NA24149_1000G_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
NA24149_1000G_PACBIO_BWAMEM_BAM_URL_PREFIX = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/PacBio_MtSinai_NIST/CSHL_bwamem_bam_GRCh37/BWA-MEM_Chr'
NA24149_1000G_PACBIO_BWAMEM_BAM_URL_SUFFIX = '_HG003_merged_11_12.sort.bam'

rule plot_pr_curve_NA24149_1000g:
    params: job_name = 'plot_pr_curve_NA24149.1000g',
            title = 'Precision Recall Curve for Reaper on NA24149: PacBio Reads vs Standard Illumina'
    input:
        reaper32_rtg = 'data/NA24149.1000g/vcfeval/reaper.pacbio.bwamem.32x.-z/{chrom}',
        illumina_rtg = 'data/NA24149.1000g/vcfeval/illumina_30x.filtered/{chrom}'
    output:
        png = 'data/plots/NA24149.1000g_prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval(['data/NA24149.1000g/vcfeval/illumina_30x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA24149.1000g/vcfeval/reaper.pacbio.bwamem.32x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 30x',
                                   'Reaper, PacBio 32x'],
                                   output.png,params.title,
                                   xlim=(0.6,1.0),ylim=(0.975,1.0))

# DOWNLOAD 60x Illumina reads
rule download_Illumina_reads_NA24149_1000g:
    params: job_name = 'DOWNLOAD_Illumina_60x_NA24149.1000g',
    output: bam = 'data/NA24149.1000g/aligned_reads/illumina/illumina.60x.bam',
    shell: 'wget {NA24149_1000G_ILLUMINA_60x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24149_1000g:
    params: job_name = 'DOWNLOAD_GIAB_BED_NA24149',
    output: 'data/NA24149.1000g/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24149_1000G_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24149_1000g:
    params: job_name = 'DOWNLOAD_GIAB_VCF_NA24149'
    output: 'data/NA24149.1000g/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24149_1000G_GIAB_VCF_URL} -O {output}'

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24149_1000g:
    params: job_name = 'download_pacbio_NA24149.1000g'
    output: bam = 'data/NA24149.1000g/aligned_reads/pacbio/pacbio.bwamem.{chrom}.32x.bam',
    shell: 'wget {NA24149_1000G_PACBIO_BWAMEM_BAM_URL_PREFIX}{wildcards.chrom}{NA24149_1000G_PACBIO_BWAMEM_BAM_URL_SUFFIX} -O {output.bam}'
