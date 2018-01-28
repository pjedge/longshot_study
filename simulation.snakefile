
import simulate_SNVs
import pysam

SIMLORD = '/home/pedge/installed/opt/python/bin/simlord'
DWGSIM = '/home/pedge/git/DWGSIM/dwgsim'
BWA = '/home/pedge/installed/bwa'
BCFTOOLS = '/opt/biotools/bcftools/bin/bcftools'
PYFAIDX = '/home/pedge/installed/opt/python/bin/faidx'
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX']

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
        fasta_hap1 = 'data/simulation/variants/ground_truth/ground_truth_hap1.fa',
        fasta_hap2 = 'data/simulation/variants/ground_truth/ground_truth_hap2.fa'
    shell:
        '''
        {BCFTOOLS} consensus -f {input.hg19} -H 1 {input.vcfgz} > {output.fasta_hap1}
        {BCFTOOLS} consensus -f {input.hg19} -H 2 {input.vcfgz} > {output.fasta_hap2}
        '''

rule join_diploid_fasta:
    params: job_name = 'generate_diploid_fasta'
    input:
        fasta_hap1 = 'data/simulation/variants/ground_truth/ground_truth_hap1.fa',
        fasta_hap2 = 'data/simulation/variants/ground_truth/ground_truth_hap2.fa'
    output:
        joined_fasta =  'data/simulation/variants/ground_truth/ground_truth.fa'
    run:
        with open(input.fasta_hap1,'r') as inf1,  open(input.fasta_hap2,'r') as inf2, open(output.joined_fasta,'w') as outf:

            printing = True
            for line in inf1:
                if line[0] == '>':
                    if line[1:].strip() in chroms:
                        printing = True
                        print(line.strip() + '_hap1', file=outf)
                    else:
                        printing = False
                else:
                    if printing:
                        print(line,end='',file=outf)

            printing = True
            for line in inf2:
                if line[0] == '>':
                    if line[1:].strip() in chroms:
                        printing = True
                        print(line.strip() + '_hap2', file=outf)
                    else:
                        printing = False
                else:
                    if printing:
                        print(line,end='',file=outf)

rule split_diploid_fasta:
    params: job_name = 'split_diploid_fasta',
            split_dir = 'data/simulation/variants/ground_truth/ground_truth_separate_chrom'
    input: joined_fasta =  'data/simulation/variants/ground_truth/ground_truth.fa'
    output: split_fasta = expand('data/simulation/variants/ground_truth/ground_truth_separate_chrom/{chrom}_hap{hap}.fa',chrom=chroms,hap=[1,2])
    shell: 'cd {params.split_dir} && {PYFAIDX} -x ../ground_truth.fa'

rule align_simulated_illumina:
    params: job_name = 'align_simulated_illumina.{cov}'
    input:
        fastq   = 'data/simulation/fastq_reads/illumina/illumina.{cov}x.fastq',
        hg19    = 'data/genomes/hg19.fa',
        hg19_ix = 'data/genomes/hg19.fa.fai',
        hg19_bwt = 'data/genomes/hg19.fa.bwt'
    output:
        bam = 'data/simulation/aligned_reads/illumina/illumina.{cov}x.bam'
    shell: '{BWA} mem -p -t 16 -T 0 {input.hg19} {input.fastq} | {SAMTOOLS} sort -@ 16 > {output.bam}'

rule align_simulated_pacbio:
    params: job_name = 'align_simulated_pacbio.{cov}'
    input:
        fastq   = 'data/simulation/fastq_reads/pacbio/pacbio.{cov}x.fastq',
        hg19    = 'data/genomes/hg19.fa',
        hg19_ix = 'data/genomes/hg19.fa.fai',
        hg19_bwt = 'data/genomes/hg19.fa.bwt'
    output:
        bam = 'data/simulation/aligned_reads/pacbio/pacbio.{cov}x.bam'
    shell: '{BWA} mem -x pacbio -t 16 -T 0 {input.hg19} {input.fastq} | {SAMTOOLS} sort -@ 16 > {output.bam}'

rule combine_reads_fastq:
    params: job_name = 'combine_{datatype}_fastq.{cov}',
    input: expand('data/simulation/fastq_reads/{{datatype}}/separate_chrom/{chrom}.hap{hap}.{{datatype}}.{{cov}}x.fastq',chrom=chroms,hap=[1,2])
    output: 'data/simulation/fastq_reads/{datatype}/{datatype}.{cov}x.fastq'
    shell: 'cat {input} > {output}'

rule simulate_illumina:
    params: job_name = 'simulate_illumina.{chrom}.hap{hap}.{cov}',
            output_prefix = 'data/simulation/fastq_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.{cov}x'
    input: diploid_fasta = 'data/simulation/variants/ground_truth/ground_truth_separate_chrom/{chrom}_hap{hap}.fa'
    output:
        fastq = 'data/simulation/fastq_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.{cov}x.fastq.gz'
    run:
        chrom_len = 0
        with pysam.FastaFile(input.diploid_fasta) as ff:
            for ref in ff.references:
                chrom_len += ff.get_reference_length(ref)

        short_read_len = 100
        haploid_cov = float(wildcards.cov) / 2.0 # the haploid coverage
        n_read_pairs = int(chrom_len * haploid_cov / (2*short_read_len))
        shell('''
        {DWGSIM} -H -z 0 -1 {short_read_len} -2 {short_read_len} -N {n_read_pairs} \
        -e 0.001 -E 0.001 -r 0 -R 0 -y 0 \
        -o 2 {input.diploid_fasta} {params.output_prefix}
        mv {params.output_prefix}.bfast.fastq.gz {params.output_prefix}.fastq.gz
        ''')

rule simulate_pacbio:
    params: job_name = 'simulate_pacbio.{chrom}.hap{hap}.{cov}',
            output_prefix = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x'
    input: diploid_fasta = 'data/simulation/variants/ground_truth/ground_truth_separate_chrom/{chrom}_hap{hap}.fa'
    output: fq = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.fastq'
    run:
        diploid_cov = int(float(wildcards.cov) / 2.0)
        shell('''
        {SIMLORD} -rr {input.diploid_fasta} \
        --coverage {diploid_cov} --no-sam {params.output_prefix}
        ''')
