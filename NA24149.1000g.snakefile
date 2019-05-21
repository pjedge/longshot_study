
NA24149_1000G_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'
NA24149_1000G_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv3.3.2/GRCh37/HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'
NA24149_1000G_PACBIO_NGMLR_BAM_URL = 'http://www.bio8.cs.hku.hk/clairvoyante/bamUsed/PacBio-HG003-hg19/all_reads.fa.giab_h003_ngmlr-0.2.3_mapped.bam'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24149_1000g:
    params: job_name = 'download_giab_bed.1000g',
    output: 'data/NA24149.1000g/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24149_1000G_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24149_1000g:
    params: job_name = 'download_giab_vcf_NA24149.1000g'
    output: 'data/NA24149.1000g/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24149_1000G_GIAB_VCF_URL} -O {output}'

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24149_1000g:
    params: job_name = 'download_pacbio_NA24149.1000g'
    output: bam = 'data/NA24149.1000g/aligned_reads/pacbio/pacbio.ngmlr.all.32x.bam',
    shell: 'wget {NA24149_1000G_PACBIO_NGMLR_BAM_URL} -O {output.bam};'
