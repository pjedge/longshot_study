
NA12878_PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
NA12878_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
NA12878_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
NA12878_Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'

rule plot_pr_curve_NA12878:
    params: job_name = 'plot_pr_curve_NA12878.1000g',
            title = 'Precision Recall Curve for Reaper on NA12878: PacBio Reads vs Standard Illumina'
    input:
        reaper30_rtg = 'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.30x.-z/{chrom}.done',
        reaper44_rtg = 'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        illumina_rtg = 'data/NA12878.1000g/vcfeval/illumina_30x.filtered/{chrom}.done'
    output:
        png = 'data/plots/NA12878.1000g_prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval(['data/NA12878.1000g/vcfeval/illumina_30x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.30x.-z/{}'.format(wildcards.chrom),
                                   'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.44x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 30x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 44x'],
                                   output.png,params.title,
                                   colors=['r','#8080ff','#3333ff'],
                                   xlim=(0.8,1.0),ylim=(0.985,1.0))

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads_NA12878:
    params: job_name = 'download_Illumina_30x_NA12878.1000g',
    output: bam = 'data/NA12878.1000g/aligned_reads/illumina/illumina.30x.bam',
    shell: 'wget {NA12878_Illumina_30x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA12878:
    params: job_name = 'download_giab_bed_NA12878.1000g',
    output: 'data/NA12878.1000g/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA12878_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA12878:
    params: job_name = 'download_giab_VCF_NA12878.1000g',
    output: 'data/NA12878.1000g/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA12878_GIAB_VCF_URL} -O {output}'

# SUBSAMPLE PACBIO BAM
rule merge_subsampled_pacbio_NA12878_1000g:
    params: job_name = 'merge_subsampled_pacbio_NA12878.1000g'
    input: bam = expand('data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.{chrom}.{{cov}}x.bam',chrom=chroms),
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.{cov,30}x.bam',
    shell: '{SAMTOOLS} merge -O bam {output.bam} {input.bam}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA12878:
    params: job_name = 'subsample_pacbio_NA12878.1000g'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.{cov}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 44.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA12878:
    params: job_name = 'split_bam_pacbio_NA12878.1000g.{chrom}'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.{chrom}.44x.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} chr{wildcards.chrom} > {output.bam}'

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA12878:
    params: job_name = 'download_pacbio_NA12878.1000g'
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    shell: 'wget {NA12878_PACBIO_BAM_URL} -O {output.bam}'
