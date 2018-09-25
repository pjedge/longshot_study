
NA12878_PACBIO_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'
NA12878_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
NA12878_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
NA12878_Illumina_30x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam'

rule plot_pr_curve_NA12878:
    params: job_name = 'plot_pr_curve_NA12878.1000g.{aligner}.{chrom}',
            title = None
    input:
        il30 = 'data/NA12878.1000g/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
        pb30 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.30x._/{chrom}',
        pb44 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.44x._/{chrom}',
        il30_cov = 'data/NA12878.1000g/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
        pb30_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage',
        pb44_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.44x.bam.median_coverage'
    output:
        png = 'data/plots/NA12878.1000g.{aligner}.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.il30, input.pb30, input.pb44],
                         ['Freebayes, Illumina {}x'.format(parse_int_file(input.il30_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb30_cov)),
                          'Longshot, PacBio {}x'.format(parse_int_file(input.pb44_cov))],
                          output.png,params.title,
                          colors=['r','#8080ff','#3333ff'],
                          xlim=(0.8,1.0),
                          ylim=(0.985,1.0))

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads_NA12878:
    params: job_name = 'download_Illumina_30x_NA12878.1000g',
    output: bam = 'data/NA12878.1000g/aligned_reads/illumina/illumina.aligned.all.30x.bam',
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

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA12878:
    params: job_name = 'split_bam_pacbio_NA12878.1000g.{chrom}'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.split_chroms/{chrom}.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} chr{wildcards.chrom} > {output.bam}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA12878:
    params: job_name = 'subsample_pacbio_NA12878.1000g'
    input: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.44x.bam',
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{chrom}.{cov,30}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 44.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA12878:
    params: job_name = 'download_pacbio_NA12878.1000g'
    output: bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
    shell: 'wget {NA12878_PACBIO_BAM_URL} -O {output.bam}'
