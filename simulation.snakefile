
import simulate_SNVs
import pysam

SIMLORD = '/home/pedge/installed/opt/python/bin/simlord'
DWGSIM = '/home/pedge/git/DWGSIM/dwgsim'
BWA = '/home/pedge/installed/bwa'
BCFTOOLS = '/opt/biotools/bcftools/bin/bcftools'
PYFAIDX = '/home/pedge/installed/opt/python/bin/faidx'
chroms = ['{}'.format(i) for i in range(1,23)] + ['X']

rule calculate_mapping_accuracy:
    params: job_name = 'calculate_mapping_accuracy{segdup_or_not}'
    input:  bam = 'data/simulation/aligned_reads/pacbio/pacbio.60x{segdup_or_not}.bam'
    output: plot = 'data/plots/chr1_simulated_60x_pacbio_mismapped_read_distribution{segdup_or_not}.png',
            acc = 'data/output/chr1_simulated_60x_pacbio_mapping_accuracy{segdup_or_not}.txt'
    run:
        acc = mapping_accuracy.mapping_accuracy(input.bam, output.plot)
        with open(output.acc,'w') as outf:
            print('mapping_accuracy={}'.format(acc),file=outf)

rule filter_chr1_segdup:
    params: job_name = 'generate_simulated_SNVs'
    input:  bam = 'data/simulation/aligned_reads/pacbio/pacbio.60x.bam',
            bed = 'genome_tracks/segmental_duplications_0.99_similar_1000g.bed'
    output: bam = 'data/simulation/aligned_reads/pacbio/pacbio.60x.segdup.bam'
    shell: '{BEDTOOLS} intersect -a {input.bam} -b {input.bed} -wa > {output.bam}'

rule plot_pr_curve_simulation:
    params: job_name = 'plot_pr_curve_simulation',
            title = 'Precision Recall Curve for Reaper on Simulated Data: PacBio Reads vs Standard Illumina'
    input:
        reaper30_rtg = 'data/simulation/vcfeval/reaper_30x.-z/{chrom}.done',
        reaper45_rtg = 'data/simulation/vcfeval/reaper_45x.-z/{chrom}.done',
        reaper60_rtg = 'data/simulation/vcfeval/reaper_60x.-z/{chrom}.done',
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
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai'
    output:
        VCF = 'data/simulation/variants/ground_truth/ground_truth.vcf'
    run:
        simulate_SNVs.simulate_SNV_VCF(input.hs37d5, output.VCF)

rule generate_diploid_fasta:
    params: job_name = 'generate_diploid_fasta'
    input:
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        vcfgz = 'data/simulation/variants/ground_truth/ground_truth.vcf.gz',
        vcf_ix = 'data/simulation/variants/ground_truth/ground_truth.vcf.gz.tbi',
    output:
        fasta_hap1 = 'data/simulation/variants/ground_truth/ground_truth_hap1.fa',
        fasta_hap2 = 'data/simulation/variants/ground_truth/ground_truth_hap2.fa'
    shell:
        '''
        {BCFTOOLS} consensus -f {input.hs37d5} -H 1 {input.vcfgz} > {output.fasta_hap1}
        {BCFTOOLS} consensus -f {input.hs37d5} -H 2 {input.vcfgz} > {output.fasta_hap2}
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
    params: job_name = 'align_simulated_illumina.{chrom}.{hap}.{cov}',
            sort_prefix = 'data/simulation/aligned_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.{cov}x.tmp'
    input:
        fastq   = 'data/simulation/fastq_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.{cov}x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = 'data/simulation/aligned_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.{cov}x.bam'
    shell: '{BWA} mem -p -t 4 -T 0 {input.hs37d5} {input.fastq} | {SAMTOOLS} sort -T {params.sort_prefix} -@ 4 > {output.bam}'

rule align_simulated_pacbio:
    params: job_name = 'align_simulated_pacbio.{chrom}.{hap}.{cov}',
            sort_prefix = 'data/simulation/aligned_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.tmp'
    input:
        fastq   = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = 'data/simulation/aligned_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.bam',
    shell: '{BWA} mem -x pacbio -t 4 -T 0 {input.hs37d5} {input.fastq} | {SAMTOOLS} sort -T {params.sort_prefix} -@ 4 > {output.bam}'

rule merge_bams:
    params: job_name = 'merge_bams.{datatype}.{cov}',
    input:  bams = expand('data/simulation/aligned_reads/{{datatype}}/separate_chrom/{chrom}.hap{hap}.{{datatype}}.{{cov}}x.bam',chrom=chroms,hap=[1,2])
    output: bam = 'data/simulation/aligned_reads/{datatype}/{datatype}.{cov}x.bam'
    shell: '{SAMTOOLS} merge -O bam {output.bam} {input.bams}'

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
