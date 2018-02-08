
rule plot_pr_curve_NA12878:
    params: job_name = 'plot_pr_curve_NA12878',
            title = 'Precision Recall Curve for Reaper on NA12878: PacBio Reads vs Standard Illumina'
    input:
        reaper30_rtg = 'data/NA12878/vcfeval/reaper_30x.-z_-C_52/{chrom}.done',
        reaper44_rtg = 'data/NA12878/vcfeval/reaper_44x.-z_-C_78/{chrom}.done',
        illumina_rtg = 'data/NA12878/vcfeval/illumina_30x.filtered/{chrom}.done'
    output:
        png = 'data/plots/NA12878_prec_recall_{chrom}.png'
    run:
        plot_vcfeval.plot_vcfeval(['data/NA12878/vcfeval/illumina_30x/filtered.{}'.format(wildcards.chrom),
                                   'data/NA12878/vcfeval/reaper_30x/{}'.format(wildcards.chrom),
                                   'data/NA12878/vcfeval/reaper_44x/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 30x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 44x'],
                                   output.png,params.title)

# DOWNLOAD 30x Illumina reads
rule download_Illumina_reads_NA12878:
    params: job_name = 'DOWNLOAD_Illumina_30x_NA12878',
    output: bam = 'data/NA12878/aligned_reads/illumina/illumina.30x.bam',
            #bai = 'data/NA12878/aligned_reads/illumina/illumina.30x.bam.bai',
    shell:
        '''
        wget {Illumina_30x_BAM_URL} -O {output.bam};
        #wget {Illumina_30x_BAI_URL} -O {output.bai};
        '''

rule extract_chrom_GIAB_VCF_NA12878:
    params: job_name = 'extract_chr_GIAB_VCF.{chrnum}',
    input:  'data/NA12878/variants/ground_truth/ground_truth.vcf'
    output: 'data/NA12878/variants/ground_truth/chr{chrnum}.vcf'
    shell: '''cat {input} | grep -P '^{chrnum}\t' > {output}'''

rule add_chr_GIAB_VCF:
    params: job_name = 'add_chr_GIAB_VCF',
    input:  gz = 'data/NA12878/variants/ground_truth/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
    output: vcf = 'data/NA12878/variants/ground_truth/ground_truth.vcf'
    shell:
        '''
        gunzip -c {input.gz} | \
        awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}' \
        > {output.vcf}
        '''

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA12878:
    params: job_name = 'DOWNLOAD_GIAB_BED_NA12878',
    output: 'data/NA12878/variants/ground_truth/region_filter.bed'
    shell: '''wget -q -O- {GIAB_HIGH_CONF_URL} | awk '{{print "chr" $0;}}' > {output}'''

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA12878:
    params: job_name = 'DOWNLOAD_GIAB_VCF_NA12878',
    output: 'data/NA12878/variants/GIAB/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
    shell: 'wget {GIAB_VCF_URL} -O {output}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA12878:
    params: job_name = 'subsample_pacbio_NA12878'
    input: bam = 'data/NA12878/aligned_reads/pacbio/pacbio.44x.bam',
    output: bam = 'data/NA12878/aligned_reads/pacbio/pacbio.{cov}x.bam',
            #bai = 'data/NA12878/aligned_reads/pacbio/pacbio.{cov}x.bam.bai'
    run:
        subsample_frac = float(wildcards.cov) / 44.0
        shell('''
        {SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam};
        #{SAMTOOLS} index {output.bam} {output.bai}
        ''')

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA12878:
    params: job_name = 'download_pacbio_NA12878'
    output: bam = 'data/NA12878/aligned_reads/pacbio/pacbio.44x.bam',
            #bai = 'data/NA12878/aligned_reads/pacbio/pacbio.44x.bam.bai'
    shell:
        '''
        wget {PACBIO_BAM_URL} -O {output.bam}
        #wget {PACBIO_BAI_URL} -O {output.bai}
        '''
