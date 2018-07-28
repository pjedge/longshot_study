
# refer to bed file with duplicated gene regions and get region for a gene
def get_gene_region(gene_bed_file, gene_name, padding=50000):
    with open(gene_bed_file, 'r') as inf:
        for line in inf:
            el = line.strip().split('\t')
            if gene_name != el[3]:
                continue

            chrom = el[0]
            start = int(el[1])+1
            end   = int(el[2])+1
            return "{}:{}-{}".format(chrom,start-padding,end+padding)

# create pacbio bamfile that is filtered for only a few low-mappability duplicated genes.
rule create_duplicated_genes_visualization_bamfile_pacbio:
    params: job_name = 'create_duplicated_genes_visualization_bamfile_pacbio.{individual}.{build}.{cov}x'
    input:  'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.CFC1.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.SMN1.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.PMS2.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.NCF1.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.STRC.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.OTOA.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.HYDIN.bam'
    output: 'data/{individual}.{build}/duplicated_gene_vis/pacbio.{cov}x.dup_genes.bam',
    shell: '{SAMTOOLS} merge -O bam {output} {input}'

# create illumina bamfile that is filtered for only a few low-mappability duplicated genes.
rule create_duplicated_genes_visualization_bamfile_illumina:
    params: job_name = 'create_duplicated_genes_visualization_bamfile_illumina.{individual}.{build}.{cov}x',
    input:  'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.CFC1.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.SMN1.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.PMS2.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.NCF1.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.STRC.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.OTOA.bam',
            'data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.HYDIN.bam'
    output: 'data/{individual}.{build}/duplicated_gene_vis/illumina.{cov}x.dup_genes.bam',
    shell: '{SAMTOOLS} merge -O bam {output} {input}'

# create pacbio bamfile that is filtered for a single duplicated gene.
rule get_single_gene_bamfile_pacbio:
    params: job_name = 'get_single_gene_bamfile_pacbio.{individual}.{build}.{cov}x',
    input:  bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.blasr.all.{cov}x.bam',
            gene_bed_file = 'genome_tracks/duplicated_genes_list_{build}.bed'
    output: bam = temp('data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/pacbio.{cov}x.{gene}.bam')
    run:
        region = get_gene_region(input.gene_bed_file, wildcards.gene)
        if wildcards.individual == 'NA12878' and wildcards.build == '1000g':
            region = 'chr' + region
        shell('{SAMTOOLS} view -hb {input.bam} {region} > {output.bam}')

# create illumina bamfile that is filtered for a single duplicated gene.
rule get_single_gene_bamfile_illumina:
    params: job_name = 'get_single_gene_bamfile_illumina.{individual}.{build}.{cov}x',
    input:  bam = 'data/{individual}.{build}/aligned_reads/illumina/illumina.{cov}x.bam',
            gene_bed_file = 'genome_tracks/duplicated_genes_list_{build}.bed'
    output: bam = temp('data/{individual}.{build}/duplicated_gene_vis/separate_genes_bams/illumina.{cov}x.{gene}.bam')
    run:
        region = get_gene_region(input.gene_bed_file, wildcards.gene)
        shell('{SAMTOOLS} view -hb {input.bam} {region} > {output.bam}')
