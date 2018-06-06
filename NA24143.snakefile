
NA24143_PACBIO_NGMLR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/all_reads.fa.giab_h004_ngmlr-0.2.3_mapped.bam'
#NA24143_PACBIO_BLASR_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg004_gr37_'
NA24143_Illumina_60x_BAM_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.60x.1.bam'
NA24143_GIAB_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh37/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'
NA24143_GIAB_HIGH_CONF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv3.3.2/GRCh37/HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed'

rule plot_pr_curve_NA24143:
    params: job_name = 'plot_pr_curve_NA24143',
            title = 'Precision Recall Curve for Reaper on NA24143: PacBio Reads vs Standard Illumina'
    input:
        reaper30_rtg = 'data/NA24143/vcfeval/reaper.pacbio.ngmlr.30x.-z/{chrom}.done',
        illumina_rtg = 'data/NA24143/vcfeval/illumina_60x.filtered/{chrom}.done'
    output:
        png = 'data/plots/NA24143_prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval(['data/NA24143/vcfeval/illumina_60x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA24143/vcfeval/reaper.pacbio.ngmlr.30x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 60x',
                                   'Reaper, PacBio 30x'],
                                   output.png,params.title,
                                   xlim=(0.6,1.0),ylim=(0.975,1.0))

# DOWNLOAD 60x Illumina reads
rule download_Illumina_reads_NA24143:
    params: job_name = 'DOWNLOAD_Illumina_60x_NA24143',
    output: bam = 'data/NA24143/aligned_reads/illumina/illumina.60x.bam',
    shell: 'wget {NA24143_Illumina_60x_BAM_URL} -O {output.bam}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_high_confidence_bed_NA24143:
    params: job_name = 'DOWNLOAD_GIAB_BED_NA24143',
    output: 'data/NA24143/variants/ground_truth/region_filter.bed'
    shell: 'wget {NA24143_GIAB_HIGH_CONF_URL} -O {output}'

# DOWNLOAD GIAB VARIANTS
rule download_GIAB_VCF_NA24143:
    params: job_name = 'DOWNLOAD_GIAB_VCF_NA24143'
    output: 'data/NA24143/variants/ground_truth/ground_truth.vcf.gz'
    shell: 'wget {NA24143_GIAB_VCF_URL} -O {output}'

# SPLIT PACBIO BAM
rule split_bam_pacbio_NA24143_NGMLR:
    params: job_name = 'split_bam_pacbio_NA24143.{chrom}'
    input: bam = 'data/NA24143/aligned_reads/pacbio/pacbio.ngmlr.all.30x.bam',
    output: bam = 'data/NA24143/aligned_reads/pacbio/pacbio.ngmlr.{chrom}.30x.bam',
    shell: '{SAMTOOLS} view -hb {input.bam} {wildcards.chrom} > {output.bam}'

# DOWNLOAD PACBIO BAM
rule download_pacbio_NA24143_NGMLR:
    params: job_name = 'download_pacbio_NA24143'
    output: bam = 'data/NA24143/aligned_reads/pacbio/pacbio.ngmlr.all.30x.bam',
    shell: 'wget {NA24143_PACBIO_NGMLR_BAM_URL} -O {output.bam}'

num_hdf5_archives = 132
hd_id = 'HG004' # NA24143 - AJ mother
hdf5_file_list = []
with open(GIAB_AJTRIO_PACBIO_HDF5_INDEX,'r') as infile:
    for line in infile:
        el = line.strip().split()
        assert(len(el) == 3)
        if el[2] != hd_id:
            continue
        hdf5_file_list.append((el[0],el[1]))

    assert(len(hdf5_file_list) == num_hdf5_archives)

# DOWNLOAD PACBIO BAM
rule bax2bam_NA24143:
    params: job_name = 'download_pacbio_NA24143'
            archive_dir = 'data/NA24143/raw_pacbio/hdf5/archive{archive_number}'
    input: 'data/NA24143/raw_pacbio/hdf5/archive{archive_number}/archive.tar.gz'
    output: 'data/NA24143/raw_pacbio/unaligned_bam/archive{archive_number}.subreads.bam'
    shell:
        '''
        cd {params.archive_dir}
        {BAX2BAM} *.bax.h5 -o ../../unaligned_bam/archive{wildcards.archive_number}
        '''

# DOWNLOAD PACBIO HDF5
rule download_pacbio_NA24143_hdf5:
    params: job_name = 'download_pacbio_NA24143_hdf5.archive{archive_number}'
            archive_dir = 'data/NA24143/raw_pacbio/hdf5/archive{archive_number}'
    output: archive = 'data/NA24143/raw_pacbio/hdf5/archive{archive_number}/archive.tar.gz',
    run:
        a = int(wildcards.archive_number)
        hdf5_url = hdf5_file_list[a][1]
        hdf5_md5 = hdf5_file_list[a][1]
        shell('''
        wget {hdf5_url} -O {output.archive}
        md5sum -c <<<"{md5_sum} *{output.archive}"
        tar -xzf {output.archive} --directory={params.archive_dir}
        ''')
