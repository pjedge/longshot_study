
NA24385_1000G_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz'
NA24385_1000G_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
NA24385_1000G_PACBIO_NGMLR_BAM_URL = 'http://www.bio8.cs.hku.hk/clairvoyante/bamUsed/PacBio-HG002-hg19/all_reads.fa.giab_h002_ngmlr-0.2.3_mapped.rg.bam'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24385_1000g:
    params: job_name = 'download_giab_bed.1000g',
    output: 'data/NA24385.1000g/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24385_1000G_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24385_1000g:
    params: job_name = 'download_giab_vcf_NA24385.1000g'
    output: 'data/NA24385.1000g/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24385_1000G_GIAB_VCF_URL} -O {output}'

# SUBSAMPLE PACBIO BAM
rule subsample_pacbio_NA24385_1000g:
    params: job_name = 'subsample_pacbio_NA24385.1000g.{aligner}.{chrom}'
    input: bam = 'data/NA24385.1000g/aligned_reads/pacbio/pacbio.{aligner}.{chrom}.69x.bam',
    output: bam = 'data/NA24385.1000g/aligned_reads/pacbio/pacbio.{aligner}.{chrom,(\d+|all)}.{cov,(20|30|40|50)}x.bam',
    run:
        subsample_frac = float(wildcards.cov) / 69.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24385_1000g:
    params: job_name = 'download_pacbio_NA24385.1000g'
    output: bam = 'data/NA24385.1000g/aligned_reads/pacbio/pacbio.ngmlr.all.69x.bam',
    shell: 'wget {NA24385_1000G_PACBIO_NGMLR_BAM_URL} -O {output.bam}'
