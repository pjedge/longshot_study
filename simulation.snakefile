
SIMLORD = '/home/pedge/installed/opt/python/bin/simlord'
WGSIM = '/opt/biotools/samtools/1.3/bin/wgsim'
BWA = '/home/pedge/installed/bwa'
BCFTOOLS = '/opt/biotools/bcftools/bin/bcftools'

HG19_LEN = 3137161264

rule plot_pr_curve_simulation:
    params: job_name = 'plot_pr_curve_simulation',
            title = 'Precision Recall Curve for Reaper on Simulated Data: PacBio Reads vs Standard Illumina'
    input:
        reaper30_rtg = 'data/simulation/vcfeval/reaper_30x.-z_-i_-B_30_-C_52/{chrom}.done',
        reaper45_rtg = 'data/simulation/vcfeval/reaper_45x.-z_-i_-B_30_-C_78/{chrom}.done',
        reaper60_rtg = 'data/simulation/vcfeval/reaper_60x.-z_-i_-B_30_-C_105/{chrom}.done',
        illumina_rtg = 'data/simulation/vcfeval/illumina_30x.filtered/{chrom}.done'
    output:
        png = 'data/plots/simulation_prec_recall_{chrom}.png'
    run:
        plot_vcfeval.plot_vcfeval(['data/simulation/vcfeval/illumina_30x/filtered.{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/reaper_30x/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/reaper_45x/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/reaper_60x/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 30x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 45x',
                                   'Reaper, PacBio 60x'],
                                   output.png,params.title)

rule generate_simulated_SNVs:
    params: job_name = 'generate_simulated_SNVs'
    input:
        hg19    = 'data/genomes/hg19.fa',
        hg19_ix = 'data/genomes/hg19.fa.fai'
    output:
        VCF = 'data/simulation/variants/ground_truth/ground_truth.vcf'
    run:
        simulate_SNVs.simulate_SNV_VCF(input.hg19, output.VCF)

rule generate_diploid_fasta:
    params: job_name = 'generate_diploid_fasta'
    input:
        hg19    = 'data/genomes/hg19.fa',
        hg19_ix = 'data/genomes/hg19.fa.fai',
        vcfgz = 'data/simulation/variants/ground_truth/ground_truth.vcf.gz',
        vcf_ix = 'data/simulation/variants/ground_truth/ground_truth.vcf.gz.tbi',
    output:
        fasta = 'data/simulation/variants/ground_truth/ground_truth.fa'
    shell:
        '''
        {BCFTOOLS} consensus -f {input.hg19} -H 1 {input.vcfgz} > {output.fasta}
        {BCFTOOLS} consensus -f {input.hg19} -H 2 {input.vcfgz} >> {output.fasta}
        '''

rule align_simulated_reads:
    params: job_name = 'align_simulated_{read_type}.{cov}'
    input:
        fastq   = 'data/simulation/fastq_reads/{read_type}/{read_type}.{cov}x.fq',
        hg19    = 'data/genomes/hg19.fa',
        hg19_ix = 'data/genomes/hg19.fa.fai'
    output:
        bam = 'data/simulation/aligned_reads/{read_type}/{read_type}.{cov}x.bam'
    run:
        pb_flag = '-x pacbio' if wildcards.read_type == 'pacbio' else ''
        shell('{BWA} mem {pb_flag} -T 0 {input.hg19} {input.fastq} | {SAMTOOLS} sort > {output.bam}')

rule simulate_reads:
    params: job_name = 'simulate_{read_type}.{cov}'
    input:
        diploid_fasta = 'data/simulation/variants/ground_truth/ground_truth.fa',
        #diploid_fasta_ix = 'data/simulation/variants/ground_truth/ground_truth.fa.fai'
    output:
        bam = 'data/simulation/fastq_reads/{read_type}/{read_type}.{cov}x.fq'
    run:

        if wildcards.read_type == 'pacbio':
            # run simlord
            diploid_cov = int(float(wildcards.cov) / 2) # divide by two since the FASTA has twice the sequences, diploid
            shell('{SIMLORD} -rr {input.diploid_fasta} --coverage {diploid_cov} --no-sam {output.bam}')
        elif wildcards.read_type == 'illumina':
            # run wgsim
            short_read_len = 100
            n_short_reads = int(HG19_LEN * float(wildcards.cov) / short_read_len)
            shell('{WGSIM} -N {n_short_reads} -e 0.001 -r 0 -R 0 -1 {short_read_len} -S 11 -d0 -e0 {input.diploid_fasta} {output.bam} /dev/null')
        else:
            print("invalid read technology")
            exit(1)
