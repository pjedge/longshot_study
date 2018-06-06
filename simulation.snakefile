
import simulate_SNVs
import pysam
import mapping_accuracy

chroms = ['{}'.format(i) for i in range(1,23)] + ['X']

# VCFeval and plotting for segmental duplications as opposed to whole genome
rule plot_pr_curve_simulation_ngmlr:
    params: job_name = 'plot_pr_curve_simulation',
            title = 'Simulated Data: Simulated PacBio Reads Aligned with NGMLR'
    input:
        reaper20_rtg = 'data/simulation/vcfeval/reaper.pacbio.ngmlr.40x.-z_-q_30/1.done',
    output:
        png = 'data/plots/simulation_prec_recall_ngmlr_1.png'
    run:
        ptf.plot_vcfeval(['data/simulation/vcfeval/reaper.pacbio.ngmlr.40x.-z_-q_30/1'],
                                   ['Reaper, PacBio 40x'],
                                   output.png,params.title,
                                   colors=['k'],
                                   xlim=(0.9,1.0),
                                   ylim=(0.9975,1.0),
                                   legendloc='lower left')

##################################################################################################################
rule plot_pr_curve_simulation_segmental_duplications_compare_mappers_no_blasr:
    params: job_name = 'plot_pr_curve_simulation_segmental_duplications_compare_mappers_no_blasr',
            title = 'Simulated 40x Data: Using Different PacBio\nMappers with Reaper in Segmental Duplications'
    input:
        bwa_q10 = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z_-q_10/{chrom}.done',
        minimap2_q10 = 'data/simulation/vcfeval_segdup/reaper.pacbio.minimap2.40x.-z_-q_10/{chrom}.done',
        ngmlr_q10 = 'data/simulation/vcfeval_segdup/reaper.pacbio.ngmlr.40x.-z_-q_10/{chrom}.done',
        bwa_q20 = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z_-q_20/{chrom}.done',
        minimap2_q20 = 'data/simulation/vcfeval_segdup/reaper.pacbio.minimap2.40x.-z_-q_20/{chrom}.done',
        ngmlr_q20 = 'data/simulation/vcfeval_segdup/reaper.pacbio.ngmlr.40x.-z_-q_20/{chrom}.done',
        bwa_q30 = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z_-q_30/{chrom}.done',
        minimap2_q30 = 'data/simulation/vcfeval_segdup/reaper.pacbio.minimap2.40x.-z_-q_30/{chrom}.done',
        ngmlr_q30 = 'data/simulation/vcfeval_segdup/reaper.pacbio.ngmlr.40x.-z_-q_30/{chrom}.done',
    output:
        png = 'data/plots/compare_mappers_no_blasr_reaper_in_segdups_simulation_{chrom,(\d+|X|Y|all)}.png'
    run:
        ptf.plot_vcfeval([input.bwa_q10[:-5], input.bwa_q20[:-5], input.bwa_q30[:-5],
                          input.minimap2_q10[:-5], input.minimap2_q20[:-5], input.minimap2_q30[:-5],
                          input.ngmlr_q10[:-5], input.ngmlr_q20[:-5], input.ngmlr_q30[:-5]],
                           ['Reaper, BWA, mapq >= 10', 'Reaper, BWA, mapq >= 20', 'Reaper, BWA, mapq >= 30',
                           'Reaper, MINIMAP2, mapq >= 10', 'Reaper, MINIMAP2, mapq >= 20', 'Reaper, MINIMAP2, mapq >= 30',
                           'Reaper, NGMLR, mapq >= 10', 'Reaper, NGMLR, mapq >= 20', 'Reaper, NGMLR, mapq >= 30',],
                           output.png,params.title,
                           colors=['#f78383','#ff4949','#ff0000',
                           '#a3a3a3','#5b5b5b','#000000',
                           '#93ff9b','#60ff6b','#00ff11'],
                           xlim=(0,1.0),
                           ylim=(0.96,1.0),
                           legendloc='lower left')

rule plot_pr_curve_simulation_segmental_duplications_compare_mappers:
    params: job_name = 'plot_pr_curve_simulation_segmental_duplications_compare_mappers',
            title = 'Simulated 40x Data: Using Different PacBio Mappers with Reaper in Segmental Duplications'
    input:
        bwa_q10 = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z_-q_10/{chrom}.done',
        blasr_q10 = 'data/simulation/vcfeval_segdup/reaper.pacbio.blasr.40x.-z_-q_10/{chrom}.done',
        minimap2_q10 = 'data/simulation/vcfeval_segdup/reaper.pacbio.minimap2.40x.-z_-q_10/{chrom}.done',
        ngmlr_q10 = 'data/simulation/vcfeval_segdup/reaper.pacbio.ngmlr.40x.-z_-q_10/{chrom}.done',
        bwa_q20 = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z_-q_20/{chrom}.done',
        blasr_q20 = 'data/simulation/vcfeval_segdup/reaper.pacbio.blasr.40x.-z_-q_20/{chrom}.done',
        minimap2_q20 = 'data/simulation/vcfeval_segdup/reaper.pacbio.minimap2.40x.-z_-q_20/{chrom}.done',
        ngmlr_q20 = 'data/simulation/vcfeval_segdup/reaper.pacbio.ngmlr.40x.-z_-q_20/{chrom}.done',
        bwa_q30 = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z_-q_30/{chrom}.done',
        blasr_q30 = 'data/simulation/vcfeval_segdup/reaper.pacbio.blasr.40x.-z_-q_30/{chrom}.done',
        minimap2_q30 = 'data/simulation/vcfeval_segdup/reaper.pacbio.minimap2.40x.-z_-q_30/{chrom}.done',
        ngmlr_q30 = 'data/simulation/vcfeval_segdup/reaper.pacbio.ngmlr.40x.-z_-q_30/{chrom}.done',
    output:
        png = 'data/plots/compare_mappers_reaper_in_segdups_simulation_{chrom,(\d+|X|Y|all)}.png'
    run:
        ptf.plot_vcfeval([input.bwa_q10[:-5], input.bwa_q20[:-5], input.bwa_q30[:-5],
                          input.blasr_q10[:-5], input.blasr_q20[:-5], input.blasr_q30[:-5],
                          input.minimap2_q10[:-5], input.minimap2_q20[:-5], input.minimap2_q30[:-5],
                          input.ngmlr_q10[:-5], input.ngmlr_q20[:-5], input.ngmlr_q30[:-5]],
                           ['Reaper, BWA, mapq >= 10', 'Reaper, BWA, mapq >= 20', 'Reaper, BWA, mapq >= 30',
                           'Reaper, BLASR, mapq >= 10', 'Reaper, BLASR, mapq >= 20', 'Reaper, BLASR, mapq >= 30',
                           'Reaper, MINIMAP2, mapq >= 10', 'Reaper, MINIMAP2, mapq >= 20', 'Reaper, MINIMAP2, mapq >= 30',
                           'Reaper, NGMLR, mapq >= 10', 'Reaper, NGMLR, mapq >= 20', 'Reaper, NGMLR, mapq >= 30',],
                           output.png,params.title,
                           colors=['#f78383','#ff4949','#ff0000',
                           '#9997fc','#6360ff','#0400ff',
                           '#a3a3a3','#5b5b5b','#000000',
                           '#93ff9b','#60ff6b','#00ff11'],
                           xlim=(0,1.0),
                           ylim=(0.98,1.0),
                           legendloc='upper right')

# VCFeval and plotting for segmental duplications as opposed to whole genome
rule plot_pr_curve_simulation_segmental_duplications:
    params: job_name = 'plot_pr_curve_simulation_segmental_duplications',
            title = 'Simulated Data: PacBio Reads vs Short Reads in Segmental Duplications'
    input:
        reaper20_rtg = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.20x.-z/{chrom}.done',
        reaper30_rtg = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.30x.-z/{chrom}.done',
        reaper40_rtg = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z/{chrom}.done',
        reaper80_rtg = 'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.80x.-z/{chrom}.done',
        illumina20_rtg = 'data/simulation/vcfeval_segdup/illumina_20x.filtered/{chrom}.done',
        illumina30_rtg = 'data/simulation/vcfeval_segdup/illumina_30x.filtered/{chrom}.done',
        illumina40_rtg = 'data/simulation/vcfeval_segdup/illumina_40x.filtered/{chrom}.done',
        illumina80_rtg = 'data/simulation/vcfeval_segdup/illumina_80x.filtered/{chrom}.done'
    output:
        png = 'data/plots/simulation_prec_recall_in_segdups_{chrom,(\d+|X|Y|all)}.png'
    run:
        ptf.plot_vcfeval(['data/simulation/vcfeval_segdup/illumina_20x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval_segdup/illumina_30x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval_segdup/illumina_40x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval_segdup/illumina_80x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.20x.-z/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.30x.-z/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.40x.-z/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval_segdup/reaper.pacbio.bwa.80x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 20x',
                                   'Freebayes, Illumina 30x',
                                   'Freebayes, Illumina 40x',
                                   'Freebayes, Illumina 80x',
                                   'Reaper, PacBio 20x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 40x',
                                   'Reaper, PacBio 80x'],
                                   output.png,params.title,
                                   colors=['#f99a9a','#fc6c6c','#ff4747','#ff0707','#9999ff','#8080ff','#6666ff','#3333ff'],
                                   xlim=(0,1.0),
                                   ylim=(0.98,1.0),
                                   legendloc='upper right')

# NOTE!!! we are filtering out indels but also MNPs which we may call as multiple SNVs
# therefore this isn't totally correct and it'd probably be better to use ROC with indels+SNVs VCF.
rule vcfeval_rtgtools_segmental_duplications:
    params: job_name = 'vcfeval_rtgtools.{dataset}.{calls_name}.{chrom}',
            region_arg = lambda wildcards: '--region={}'.format(wildcards.chrom) if wildcards.chrom != 'all' else ''
    input:  calls_vcf = 'data/{dataset}/variants/{calls_name}/{chrom}.vcf.gz',
            calls_ix = 'data/{dataset}/variants/{calls_name}/{chrom}.vcf.gz.tbi',
            ground_truth = 'data/{dataset}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz',
            ground_truth_ix = 'data/{dataset}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz.tbi',
            region_filter ='genome_tracks/segmental_duplications_0.99_similar_1000g.bed',
            tg_sdf = 'data/genomes/1000g_v37_phase2.sdf'
    output: done = 'data/{dataset}/vcfeval_segdup/{calls_name}/{chrom,(\d+|X|Y|all)}.done'
    shell:
        '''
        rm -rf data/{wildcards.dataset}/vcfeval_segdup/{wildcards.calls_name}/{wildcards.chrom}
        {RTGTOOLS} RTG_MEM=12g vcfeval \
        {params.region_arg} \
        -c {input.calls_vcf} \
        -b {input.ground_truth} \
        -e {input.region_filter} \
        -t {input.tg_sdf} \
        -o data/{wildcards.dataset}/vcfeval_segdup/{wildcards.calls_name}/{wildcards.chrom};
        cp data/{wildcards.dataset}/vcfeval_segdup/{wildcards.calls_name}/{wildcards.chrom}/done {output.done};
        '''
##################################################################################################################

rule calculate_mapping_accuracy:
    params: job_name = 'calculate_mapping_accuracy{segdup_or_not}{chrom}'
    input:  bam = 'data/simulation/aligned_reads/pacbio/pacbio.bwa.{chrom}.60x{segdup_or_not}bam'
    output: plot = 'data/plots/chr{chrom}_simulated_60x_pacbio_mismapped_read_distribution{segdup_or_not}png',
            acc = 'data/output/chr{chrom}_simulated_60x_pacbio_mapping_accuracy{segdup_or_not}txt'
    run:
        acc = mapping_accuracy.mapping_accuracy(input.bam, output.plot)
        with open(output.acc,'w') as outf:
            print('mapping_accuracy={}'.format(acc),file=outf)

rule filter_chr1_segdup:
    params: job_name = 'filter_chr1_segdup.{chrom}'
    input:  bam = 'data/simulation/aligned_reads/pacbio/pacbio.bwa.{chrom}.60x.bam',
            bed = 'genome_tracks/segmental_duplications_0.99_similar_1000g.bed'
    output: bam = 'data/simulation/aligned_reads/pacbio/pacbio.bwa.{chrom}.60x.segdup.bam'
    shell: '{BEDTOOLS} intersect -a {input.bam} -b {input.bed} -wa > {output.bam}'

rule plot_pr_curve_simulation:
    params: job_name = 'plot_pr_curve_simulation',
            title = 'Precision Recall Curve for Reaper on Simulated Data: PacBio Reads vs Standard Illumina'
    input:
        reaper20_rtg = 'data/simulation/vcfeval/reaper.pacbio.bwa.20x.-z/{chrom}.done',
        reaper30_rtg = 'data/simulation/vcfeval/reaper.pacbio.bwa.30x.-z/{chrom}.done',
        reaper40_rtg = 'data/simulation/vcfeval/reaper.pacbio.bwa.40x.-z/{chrom}.done',
        reaper80_rtg = 'data/simulation/vcfeval/reaper.pacbio.bwa.80x.-z/{chrom}.done',
        illumina20_rtg = 'data/simulation/vcfeval/illumina_20x.filtered/{chrom}.done',
        illumina30_rtg = 'data/simulation/vcfeval/illumina_30x.filtered/{chrom}.done',
        illumina40_rtg = 'data/simulation/vcfeval/illumina_40x.filtered/{chrom}.done',
        illumina80_rtg = 'data/simulation/vcfeval/illumina_80x.filtered/{chrom}.done'
    output:
        png = 'data/plots/simulation_prec_recall_{chrom,(\d+|X|Y|all)}.png'
    run:
        ptf.plot_vcfeval(['data/simulation/vcfeval/illumina_20x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/illumina_30x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/illumina_40x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/illumina_80x.filtered/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/reaper.pacbio.bwa.20x.-z/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/reaper.pacbio.bwa.30x.-z/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/reaper.pacbio.bwa.40x.-z/{}'.format(wildcards.chrom),
                                   'data/simulation/vcfeval/reaper.pacbio.bwa.80x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 20x',
                                   'Freebayes, Illumina 30x',
                                   'Freebayes, Illumina 40x',
                                   'Freebayes, Illumina 80x',
                                   'Reaper, PacBio 20x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 40x',
                                   'Reaper, PacBio 80x'],
                                   output.png,params.title,
                                   colors=['#f99a9a','#fc6c6c','#ff4747','#ff0707','#9999ff','#8080ff','#6666ff','#3333ff'],
                                   xlim=(0.8,1.0),
                                   ylim=(0.99,1.0))

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
        bam = 'data/simulation/aligned_reads/illumina/separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.illumina.{cov,\d+}x.bam'
    shell: '{BWA} mem -p -t 4 -T 0 {input.hs37d5} {input.fastq} | {SAMTOOLS} sort -T {params.sort_prefix} -@ 4 > {output.bam}'

rule align_simulated_pacbio_bwa:
    params: job_name = 'align_simulated_pacbio_bwa.{chrom}.{hap}.{cov}',
            sort_prefix = 'data/simulation/aligned_reads/pacbio/bwa_separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.tmp'
    input:
        fastq   = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = 'data/simulation/aligned_reads/pacbio/bwa_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.{cov,\d+}x.bam',
    shell: '{BWA} mem -x pacbio -t 4 -T 0 {input.hs37d5} {input.fastq} | {SAMTOOLS} sort -T {params.sort_prefix} -@ 4 > {output.bam}'

rule convert_simulated_fastqs_pacbio_format:
    params: job_name = 'convert_simulated_fastqs_pacbio_format.chr{chrom}.hap{hap}.{cov}x'
    input:  'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.fastq'
    output: 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.pacbio_format.fa'
    run:
        with pysam.FastxFile(input[0]) as infile, open(output[0],'w') as outfile:
            for i,record in enumerate(infile):
                print(">{}/{}/0_{}".format(record.name,i,len(record.sequence)),file=outfile)
                print(record.sequence, file=outfile)

rule align_simulated_pacbio_blasr:
    params: job_name = 'align_simulated_pacbio_blasr.{chrom}.{hap}.{cov}'
    input:
        fastq   = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.pacbio_format.fa',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_sa = 'data/genomes/hs37d5.fa.sawriter.sa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output: sam = 'data/simulation/aligned_reads/pacbio/blasr_separate_chrom_unsorted/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.{cov,\d+}x.bam',
    shell: '{BLASR} {input.fastq} {input.hs37d5} --sa {input.hs37d5_sa} --nproc 4 --bam --out {output}'

rule sort_simulated_pacbio_blasr:
    params: job_name = 'sort_simulated_pacbio_blasr.{chrom}.{hap}.{cov}',
            sort_prefix = 'data/simulation/aligned_reads/pacbio/blasr_separate_chrom_unsorted/{chrom}.hap{hap}.pacbio.{cov}x.tmp'
    input: 'data/simulation/aligned_reads/pacbio/blasr_separate_chrom_unsorted/{chrom}.hap{hap}.pacbio.{cov}x.bam'
    output: 'data/simulation/aligned_reads/pacbio/blasr_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.{cov,\d+}x.bam',
    shell: '{SAMTOOLS} sort -T {params.sort_prefix} -@ 4 {input} > {output}'

rule make_blasr_suffix_array:
    params: job_name = 'make_blasr_suffix_array.{genome}'
    input: 'data/genomes/{genome}.fa'
    output: 'data/genomes/{genome}.fa.sawriter.sa'
    shell: '{SAWRITER} {output} {input}'

rule align_simulated_pacbio_minimap2:
    params: job_name = 'align_simulated_pacbio_minimap2.{chrom}.{hap}.{cov}',
            sort_prefix = 'data/simulation/aligned_reads/pacbio/minimap2_separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.tmp'
    input:
        fastq   = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_mmi    = 'data/genomes/hs37d5.fa.mmi',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = 'data/simulation/aligned_reads/pacbio/minimap2_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.{cov,\d+}x.bam',
    shell: '{MINIMAP2} -ax map-pb {input.hs37d5_mmi} {input.fastq} | {SAMTOOLS} view -hb | {SAMTOOLS} sort -T {params.sort_prefix} -@ 4 > {output.bam}'

rule align_simulated_pacbio_ngmlr:
    params: job_name = 'align_simulated_pacbio_ngmlr.{chrom}.{hap}.{cov}',
            sort_prefix = 'data/simulation/aligned_reads/pacbio/ngmlr_separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.tmp'
    input:
        fastq   = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom}.hap{hap}.pacbio.{cov}x.fastq',
        hs37d5    = 'data/genomes/hs37d5.fa',
        hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
        hs37d5_bwt = 'data/genomes/hs37d5.fa.bwt'
    output:
        bam = 'data/simulation/aligned_reads/pacbio/ngmlr_separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.{cov,\d+}x.bam',
    shell: '{NGMLR} -t 4 -x pacbio -r {input.hs37d5} -q {input.fastq} | {SAMTOOLS} view -hb | {SAMTOOLS} sort -T {params.sort_prefix} -@ 4 > {output.bam}'

rule merge_illumina_chrom_bams:
    params: job_name = 'merge_illumina_chrom_bams.{cov}',
    input: expand('data/simulation/aligned_reads/illumina/illumina.bwa.{chrom}.{{cov}}x.bam', chrom=chroms)
    output: 'data/simulation/aligned_reads/illumina/illumina.{cov,\d+}x.bam'
    shell: '{SAMTOOLS} merge -O bam {output} {input}'

################################################################################
# IMPORTANT
# we are currently simulating pacbio reads from a single chromosome, mapping those to
# the genome, and then calling variants on that single chromosome

# this saves time having to simulate and map the entire genome, but it makes
# the results less valid because we are excluding the effect of reads generated
# from other chromosomes that mismap into the chromosome we're calling variants in.
################################################################################

rule merge_haplotype_bams:
    params: job_name = 'merge_haplotype_bams.{datatype}.{chrom}.{cov}.{aligner}',
    input: 'data/simulation/aligned_reads/{datatype}/{aligner}_separate_chrom/{chrom}.hap1.{datatype}.{cov}x.bam',
           'data/simulation/aligned_reads/{datatype}/{aligner}_separate_chrom/{chrom}.hap2.{datatype}.{cov}x.bam'
    output: 'data/simulation/aligned_reads/{datatype}/{datatype}.{aligner}.{chrom,(\d+|X|Y|all)}.{cov,\d+}x.bam'
    shell: '{SAMTOOLS} merge -O bam {output} {input}'

rule simulate_illumina:
    params: job_name = 'simulate_illumina.{chrom}.hap{hap}.{cov}',
            output_prefix = 'data/simulation/fastq_reads/illumina/separate_chrom/{chrom}.hap{hap}.illumina.{cov}x'
    input: diploid_fasta = 'data/simulation/variants/ground_truth/ground_truth_separate_chrom/{chrom}_hap{hap}.fa'
    output:
        fastq = 'data/simulation/fastq_reads/illumina/separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap,(1|2)}.illumina.{cov,\d+}x.fastq.gz'
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
    output: fq = 'data/simulation/fastq_reads/pacbio/separate_chrom/{chrom,(\d+|X|Y|all)}.hap{hap}.pacbio.{cov,\d+}x.fastq'
    run:
        diploid_cov = int(float(wildcards.cov) / 2.0)
        shell('''
        {SIMLORD} -rr {input.diploid_fasta} \
        --coverage {diploid_cov} --no-sam {params.output_prefix}
        ''')

rule copy_simulated_region_filter:
    params: job_name = 'copy_simulated_region_filter'
    input: 'genome_tracks/whole_genome_1000g.bed',
    output: 'data/simulation/variants/ground_truth/region_filter.bed',
    shell: 'cp {input} {output}'
