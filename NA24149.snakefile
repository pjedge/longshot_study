
NA24149_PACBIO_NGMLR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/all_reads.fa.giab_h003_ngmlr-0.2.3_mapped.bam'
#NA24149_PACBIO_BLASR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg003_gr37_'
NA24149_Illumina_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.60x.1.bam'
NA24149_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'
NA24149_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'

rule plot_pr_curve_NA24149:
    params: job_name = 'plot_pr_curve_NA24149',
            title = 'Precision Recall Curve for Reaper on NA24149: PacBio Reads vs Standard Illumina'
    input:
        reaper32_rtg = 'data/NA24149/vcfeval/reaper.pacbio.ngmlr.32x.-z/{chrom}.done',
        illumina_rtg = 'data/NA24149/vcfeval/illumina_60x.filtered/{chrom}.done'
    output:
        png = 'data/plots/NA24149_prec_recall_{chrom}.png'
    run:
        plot_vcfeval.plot_vcfeval(['data/NA24149/vcfeval/illumina_60x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA24149/vcfeval/reaper.pacbio.ngmlr.32x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 60x',
                                   'Reaper, PacBio 32x'],
                                   output.png,params.title)

# DOWNLOAD 60x Illumina reads
rule download_Illumina_reads_NA24149:
    params: job_name = 'DOWNLOAD_Illumina_60x_NA24149',
    output: bam = 'data/NA24149/aligned_reads/illumina/illumina.60x.bam',
    shell: 'wget {NA24149_Illumina_60x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24149:
    params: job_name = 'DOWNLOAD_GIAB_BED_NA24149',
    output: 'data/NA24149/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24149_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24149:
    params: job_name = 'DOWNLOAD_GIAB_VCF_NA24149'
    output: 'data/NA24149/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24149_GIAB_VCF_URL} -O {output}'

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA24149_NGMLR:
    params: job_name = 'split_bam_pacbio_NA24149.{chrom}'
    input: bam = 'data/NA24149/aligned_reads/pacbio/pacbio.ngmlr.all.32x.bam',
    output: bam = 'data/NA24149/aligned_reads/pacbio/pacbio.ngmlr.{chrom}.32x.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} chr{wildcards.chrom} > {output.bam}'

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24149_NGMLR:
    params: job_name = 'download_pacbio_NA24149'
    output: bam = 'data/NA24149/aligned_reads/pacbio/pacbio.ngmlr.all.32x.bam',
    shell: 'wget {NA24149_PACBIO_NGMLR_BAM_URL} -O {output.bam}'

# DOWNLOAD PACBIO BAM
#rule download_pacbio_NA24149:
#    params: job_name = 'download_pacbio_NA24149.{chrom}'
#    output: bam = 'data/NA24149/aligned_reads/pacbio/pacbio.blasr.{chrom}.32x.bam',
#    shell: 'wget {NA24149_PACBIO_BLASR_BAM_URL}{wildcards.chrom}.bam -O {output.bam}'
