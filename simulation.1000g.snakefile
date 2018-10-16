
import simulate_SNVs
import pysam
import mapping_accuracy

rule vcfeval_rtgtools_segmental_duplications:
    params: job_name = 'vcfeval_rtgtools_segdup.{dataset}.1000g.{calls_name}.{chrom}',
            region_arg = lambda wildcards: '--region={}'.format(wildcards.chrom) if wildcards.chrom != 'all' else ''
    input:  calls_vcf = 'data/{dataset}.1000g/variants/{calls_name}/{chrom}.vcf.gz',
            calls_ix = 'data/{dataset}.1000g/variants/{calls_name}/{chrom}.vcf.gz.tbi',
            ground_truth = 'data/{dataset}.1000g/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz',
            ground_truth_ix = 'data/{dataset}.1000g/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz.tbi',
            region_filter ='genome_tracks/segmental_duplications_0.95_similar_1000g.bed',
            tg_sdf = 'data/genomes/1000g.sdf'
    output: dir = directory('data/{dataset}.1000g/vcfeval_segdup/{calls_name}/{chrom,(\d+|X|Y|all)}')
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcfeval \
        {params.region_arg} \
        -c {input.calls_vcf} \
        -b {input.ground_truth} \
        -e {input.region_filter} \
        -t {input.tg_sdf} \
        -o {output.dir}
        '''
##################################################################################################################

rule calculate_mapping_accuracy:
    params: job_name = 'calculate_mapping_accuracy{segdup_or_not}{chrom}'
    input:  bam = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.blasr.{chrom}.60x{segdup_or_not}bam'
    output: plot = 'data/plots/chr{chrom}_simulated_60x_pacbio_mismapped_read_distribution{segdup_or_not}png',
            acc = 'data/output/chr{chrom}_simulated_60x_pacbio_mapping_accuracy{segdup_or_not}txt'
    run:
        acc = mapping_accuracy.mapping_accuracy(input.bam, output.plot)
        with open(output.acc,'w') as outf:
            print('mapping_accuracy={}'.format(acc),file=outf)

rule combine_segdup_stats:
    params: job_name = 'combine_segdup_stats'
    input:  expand('data/simulation.1000g/aligned_reads/pacbio/segdup95_stats/{chrom}.stats.bed',chrom=chroms)
    output: 'data/simulation.1000g/aligned_reads/pacbio/segdup95_stats/all.stats.bed'
    shell: 'cat {input} > {output}'

rule get_segdup_stats:
    params: job_name = 'get_segdup_stats.{chrom}'
    input:  bam = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.blasr.{chrom}.60x.filtered.mapq30.segdup95.bam',
            bai = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.blasr.{chrom}.60x.filtered.mapq30.segdup95.bam.bai',
            bed = 'genome_tracks/segmental_duplications_0.95_similar_1000g.bed'
    output: bed = 'data/simulation.1000g/aligned_reads/pacbio/segdup95_stats/{chrom}.stats.bed'
    run:
        mapping_accuracy.mapping_accuracy_and_completeness_segdups(bamfile=input.bam,
                                                                   bedfile=input.bed,
                                                                   outfile=output.bed,
                                                                   chrom_filter=wildcards.chrom,
                                                                   min_mapq=30,
                                                                   delta=5000)

rule filter_simulation_segdup_blasr:
    params: job_name = 'filter_simulation_segdup.blasr.{chrom}'
    input:  bam = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.blasr.{chrom}.60x.bam',
            bed = 'genome_tracks/segmental_duplications_0.95_similar_1000g.bed'
    output: bam = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.blasr.{chrom}.60x.filtered.mapq30.segdup95.bam'
    shell: '{SAMTOOLS} view -hb -F 3844 -q 30 {input.bam} -L {input.bed} > {output.bam}'

####################################################################################################################

rule generate_simulated_SNVs:
    params: job_name = 'generate_simulated_SNVs'
    input:
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai'
    output:
        VCF = 'data/simulation.1000g/variants/ground_truth/ground_truth.vcf'
    run:
        simulate_SNVs.simulate_SNV_VCF(input.hs37d5, output.VCF)

rule generate_diploid_fasta:
    params: job_name = 'generate_diploid_fasta'
    input:
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        vcfgz = 'data/simulation.1000g/variants/ground_truth/ground_truth.vcf.gz',
        vcf_ix = 'data/simulation.1000g/variants/ground_truth/ground_truth.vcf.gz.tbi',
    output:
        fasta_hap1 = 'data/simulation.1000g/variants/ground_truth/ground_truth_hap1.fa',
        fasta_hap2 = 'data/simulation.1000g/variants/ground_truth/ground_truth_hap2.fa'
    shell:
        '''
        {BCFTOOLS} consensus -f {input.hs37d5} -H 1 {input.vcfgz} > {output.fasta_hap1}
        {BCFTOOLS} consensus -f {input.hs37d5} -H 2 {input.vcfgz} > {output.fasta_hap2}
        '''

rule join_diploid_fasta:
    params: job_name = 'generate_diploid_fasta'
    input:
        fasta_hap1 = 'data/simulation.1000g/variants/ground_truth/ground_truth_hap1.fa',
        fasta_hap2 = 'data/simulation.1000g/variants/ground_truth/ground_truth_hap2.fa'
    output:
        joined_fasta =  'data/simulation.1000g/variants/ground_truth/ground_truth.fa'
    run:
        with open(input.fasta_hap1,'r') as inf1,  open(input.fasta_hap2,'r') as inf2, open(output.joined_fasta,'w') as outf:

            printing = True
            for line in inf1:
                if line[0] == '>':
                    chrname = line[1:].strip().split(' ')[0]
                    if chrname in chroms:
                        printing = True
                        print('>' + chrname + '_hap1', file=outf)
                    else:
                        printing = False
                else:
                    if printing:
                        print(line,end='',file=outf)

            printing = True
            for line in inf2:
                if line[0] == '>':
                    chrname = line[1:].strip().split(' ')[0]
                    if chrname in chroms:
                        printing = True
                        print('>' + chrname + '_hap2', file=outf)
                    else:
                        printing = False
                else:
                    if printing:
                        print(line,end='',file=outf)

rule split_diploid_fasta:
    params: job_name = 'split_diploid_fasta',
            split_dir = 'data/simulation.1000g/variants/ground_truth/ground_truth_separate_chrom'
    input: joined_fasta =  'data/simulation.1000g/variants/ground_truth/ground_truth.fa'
    output: split_fasta = expand('data/simulation.1000g/variants/ground_truth/ground_truth_separate_chrom/{chrom}_hap{hap}.fa',chrom=chroms,hap=[1,2])
    shell: 'cd {params.split_dir} && {PYFAIDX} -x ../ground_truth.fa'

rule align_simulated_illumina:
    params: job_name = 'align_simulated_illumina.{chrom}.{hap}.60x',
            sort_prefix = 'data/simulation.1000g/aligned_reads/illumina/bwa_separate_chrom/{chrom}.hap{hap}.illumina.60x.tmp'
    input:
        fastq   = 'data/simulation.1000g/fastq_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.60x.fastq.gz',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = temp('data/simulation.1000g/aligned_reads/illumina/aligned_separate_chrom/{chrom,(\d+|all)}.hap{hap}.illumina.60x.bam')
    shell: '{BWA} mem -p -t 4 -T 0 {input.hs37d5} {input.fastq} | {SAMTOOLS} sort -T {params.sort_prefix} -@ 4 > {output.bam}'

rule align_simulated_pacbio_bwa:
    params: job_name = 'align_simulated_pacbio_bwa.{chrom}.{hap}.60x',
            sort_prefix = 'data/simulation.1000g/aligned_reads/pacbio/bwa_separate_chrom/{chrom}.hap{hap}.pacbio.60x.tmp'
    input:
        fastq   = 'data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.60x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = temp('data/simulation.1000g/aligned_reads/pacbio/bwa_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.60x.bam'),
    shell: '{BWA} mem -x pacbio -t 16 -T 0 {input.hs37d5} {input.fastq} | {SAMTOOLS} sort -T {params.sort_prefix} -@ 16 > {output.bam}'

rule convert_simulated_fastqs_pacbio_format:
    params: job_name = 'convert_simulated_fastqs_pacbio_format.chr{chrom}.hap{hap}.60x'
    input:  'data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.60x.fastq'
    output: temp('data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.60x.pacbio_format.fa')
    run:
        with pysam.FastxFile(input[0]) as infile, open(output[0],'w') as outfile:
            for i,record in enumerate(infile):
                print(">{}/{}/0_{}".format(record.name,i,len(record.sequence)),file=outfile)
                print(record.sequence, file=outfile)

rule align_simulated_pacbio_blasr:
    params: job_name = 'align_simulated_pacbio_blasr.{chrom}.{hap}.60x'
    input:
        fastq   = 'data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.60x.pacbio_format.fa',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_sa = 'data/genomes/hs37d5.fa.sawriter.sa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output: bam = temp('data/simulation.1000g/aligned_reads/pacbio/blasr_separate_chrom_unsorted/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.60x.bam'),
    shell: '{BLASR} {input.fastq} {input.hs37d5} --sa {input.hs37d5_sa} --nproc 16 --bestn 1 --bam --out {output}'

rule sort_simulated_pacbio_blasr:
    params: job_name = 'sort_simulated_pacbio_blasr.{chrom}.{hap}.60x',
            sort_prefix = 'data/simulation.1000g/aligned_reads/pacbio/blasr_separate_chrom_unsorted/{chrom}.hap{hap}.pacbio.60x.tmp'
    input: 'data/simulation.1000g/aligned_reads/pacbio/blasr_separate_chrom_unsorted/{chrom}.hap{hap}.pacbio.60x.bam'
    output: temp('data/simulation.1000g/aligned_reads/pacbio/blasr_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.60x.bam'),
    shell: '{SAMTOOLS} sort -T {params.sort_prefix} -@ 16 {input} > {output}'

rule make_blasr_suffix_array:
    params: job_name = 'make_blasr_suffix_array.{genome}'
    input: 'data/genomes/{genome}.fa'
    output: 'data/genomes/{genome}.fa.sawriter.sa'
    shell: '{SAWRITER} {output} {input}'

rule align_simulated_pacbio_minimap2:
    params: job_name = 'align_simulated_pacbio_minimap2.{chrom}.{hap}.60x',
            sort_prefix = 'data/simulation.1000g/aligned_reads/pacbio/minimap2_separate_chrom/{chrom}.hap{hap}.pacbio.60x.tmp'
    input:
        fastq   = 'data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.60x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_mmi    = 'data/genomes/hs37d5.fa.mmi',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = temp('data/simulation.1000g/aligned_reads/pacbio/minimap2_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.60x.bam'),
    shell: '{MINIMAP2} -t 16 -ax map-pb {input.hs37d5_mmi} {input.fastq} | {SAMTOOLS} view -hb | {SAMTOOLS} sort -T {params.sort_prefix} -@ 16 > {output.bam}'

rule align_simulated_pacbio_ngmlr:
    params: job_name = 'align_simulated_pacbio_ngmlr.{chrom}.{hap}.60x',
            sort_prefix = 'data/simulation.1000g/aligned_reads/pacbio/ngmlr_separate_chrom/{chrom}.hap{hap}.pacbio.60x.tmp'
    input:
        fastq   = 'data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.60x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = temp('data/simulation.1000g/aligned_reads/pacbio/ngmlr_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.60x.bam'),
    shell: '{NGMLR} -t 16 -x pacbio -r {input.hs37d5} -q {input.fastq} | {SAMTOOLS} view -hb | {SAMTOOLS} sort -T {params.sort_prefix} -@ 16 > {output.bam}'

# SUBSAMPLE PACBIO BAM
rule subsample_simulated_reads:
    params: job_name = 'subsample_simulated_{tech}.{aligner}.{cov}x'
    input: bam = 'data/simulation.1000g/aligned_reads/{tech}/{tech}.{aligner}.all.60x.bam',
    output: bam = 'data/simulation.1000g/aligned_reads/{tech}/{tech}.{aligner}.all.{cov,(20|30|40)}x.bam'
    run:
        subsample_frac = float(wildcards.cov) / 60.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

# SPLIT PACBIO BAM
rule split_simulated_pacbio_bam:
    params: job_name = 'split_simulated_pacbio_bam.{cov}.{chrom}'
    input: bam = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam',
           bai = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam.bai'
    output: bam = 'data/simulation.1000g/aligned_reads/pacbio/pacbio.{aligner}.{chrom,(\d+)}.{cov,(20|30|40|60)}x.bam'
    shell: '{SAMTOOLS} view -hb {input.bam} {wildcards.chrom} > {output.bam}'

rule merge_simulated_bams:
    params: job_name = 'merge_simulated_bams.{datatype}.60x.{aligner}',
    input: expand('data/simulation.1000g/aligned_reads/{{datatype}}/{{aligner}}_separate_chrom/{chrom}.hap{h}.{{datatype}}.60x.bam',h=[1,2],chrom=chroms)
    output: 'data/simulation.1000g/aligned_reads/{datatype}/{datatype}.{aligner}.all.60x.bam'
    shell: '{SAMTOOLS} merge -O bam {output} {input}'

rule simulate_illumina:
    params: job_name = 'simulate_illumina.{chrom}.hap{hap}.60x',
            output_prefix = 'data/simulation.1000g/fastq_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.60x'
    input: diploid_fasta = 'data/simulation.1000g/variants/ground_truth/ground_truth_separate_chrom/{chrom}_hap{hap}.fa'
    output:
        fastq = 'data/simulation.1000g/fastq_reads/illumina/separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap,(1|2)}.illumina.60x.fastq.gz'
    run:
        chrom_len = 0
        with pysam.FastaFile(input.diploid_fasta) as ff:
            for ref in ff.references:
                seq = ff.fetch(ref)
                chrom_len = 0
                for char in seq:
                    if char != 'N' and char != 'n':
                        chrom_len += 1
                #chrom_len += ff.get_reference_length(ref)

        short_read_len = 100
        haploid_cov = 30 #float(wildcards.cov) / 2.0 # the haploid coverage
        n_read_pairs = int(chrom_len * haploid_cov / (2*short_read_len))
        shell('''
        {DWGSIM} -H -z 0 -1 {short_read_len} -2 {short_read_len} -N {n_read_pairs} \
        -e 0.001 -E 0.001 -r 0 -R 0 -y 0 \
        -o 2 {input.diploid_fasta} {params.output_prefix}
        mv {params.output_prefix}.bfast.fastq.gz {params.output_prefix}.fastq.gz
        ''')

rule simulate_pacbio:
    params: job_name = 'simulate_pacbio.{chrom}.hap{hap}.60x',
            output_prefix = 'data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.60x'
    input: diploid_fasta = 'data/simulation.1000g/variants/ground_truth/ground_truth_separate_chrom/{chrom}_hap{hap}.fa'
    output: fq = 'data/simulation.1000g/fastq_reads/pacbio/separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.60x.fastq'
    run:
        #this is the coverage you need to generate reads at for a diploid fasta,
        # in order to get 60x when you align to haploid reference
        diploid_cov = 30 #int(float(wildcards.cov) / 2.0)
        shell('''
        {SIMLORD} -mp 1 -rr {input.diploid_fasta} \
        --coverage {diploid_cov} --no-sam {params.output_prefix}
        ''')

rule copy_simulated_region_filter:
    params: job_name = 'copy_simulated_region_filter'
    input: 'genome_tracks/whole_genome_1000g.bed',
    output: 'data/simulation.1000g/variants/ground_truth/region_filter.bed',
    shell: 'cp {input} {output}'
