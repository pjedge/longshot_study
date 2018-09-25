
import re

DUP_GENES_FOR_TABLE = [
'NOTCH2',
'CFC1',
'OCLN',
'RPS17',
'SMN1',
'CYP21A2',
'NCF1',
'PMS2',
'ADAMTSL2',
'STRC',
'ABCC6',
'HYDIN',
'OTOA',
'PKD1',
#'IKBKG',
#'OPN1MW'
]

# takes a gene name and path to a gencode bed file
# returns the region (chrom,start,end) (0-based,exclusive) of the gene
def search_gencode_bed(gene_name, table_file):

    region = None
    gene_name_re = re.compile("gene_name=([^;]+);")
    with open(table_file,'r') as infile:
        for line in infile:
            el = line.strip().split('\t')
            chrom = el[0]
            start = int(el[1])
            end   = int(el[2])
            data  = el[9]
            gene = gene_name_re.findall(data)[0]

            if gene == gene_name:
                assert(region == None)
                region = (chrom, start, end)

    assert(region != None)

    return region

# takes a gene name and path to a refGene table file
# returns the region (chrom,start,end) (0-based,exclusive) of the gene
#def search_refgene_table(gene_name, table_file):
#
#    chrom = None
#    start = None
#    end = None
#
#    with open(table_file,'r') as infile:
#        for line in infile:
#            el = line.strip().split('\t')
#            name = el[12]
#            if name == gene_name:
#                #if region != None:
#                #    print(gene_name," has an extra entry")
#                #    continue
#                #assert(region == None)
#                entry_chrom = el[2]
#                entry_start = int(el[4])
#                entry_end = int(el[5])
#
#                if chrom == None:
#                    chrom = entry_chrom
#
#                if chrom != entry_chrom:
#                    import pdb; pdb.set_trace()
#                assert(chrom == entry_chrom)
#
#                if start == None or entry_start < start:
#                    start = entry_start
#
#                if end == None or entry_end > end:
#                    end = entry_end

    #if region == None:
    #    import pdb; pdb.set_trace()

#    assert(chrom != None)
#    assert(start != None)
#    assert(end != None)

#    return (chrom, start, end)

def parse_region_file(region_file):
    with open(region_file,'r') as inf:
        chrom, start_end = inf.read().strip().split(":")
        start,end = [int(x) for x in start_end.split("-")]
        return (chrom,start,end)

rule make_gene_table:
    params: job_name = 'make_gene_table.{id}.{build}.pb{pcov}x'
    input:  regions = expand('data/{{id}}.{{build}}/gene_mappability_pb{{pcov}}x/{gene_name}/region.txt',gene_name=DUP_GENES_FOR_TABLE),
            pb_SNV_stats = expand('data/{{id}}.{{build}}/gene_mappability_pb{{pcov}}x/{gene_name}/pacbio_SNVs.vcf.stats',gene_name=DUP_GENES_FOR_TABLE),
            il_SNV_stats = expand('data/{{id}}.{{build}}/gene_mappability_pb{{pcov}}x/{gene_name}/illumina_SNVs.vcf.stats',gene_name=DUP_GENES_FOR_TABLE),
            shared_SNV_stats = expand('data/{{id}}.{{build}}/gene_mappability_pb{{pcov}}x/{gene_name}/illumina_intersect_pacbio_SNVs.vcf.stats',gene_name=DUP_GENES_FOR_TABLE),
            pb_frac_map = expand('data/{{id}}.{{build}}/gene_mappability_pb{{pcov}}x/{gene_name}/pacbio_frac_mappable.txt',gene_name=DUP_GENES_FOR_TABLE),
            il_frac_map = expand('data/{{id}}.{{build}}/gene_mappability_pb{{pcov}}x/{gene_name}/illumina_frac_mappable.txt',gene_name=DUP_GENES_FOR_TABLE),
    output: tex = 'data/output/duplicated_genes_table.{id}.{build}.pb{pcov}x.tex'
    run:
        ptf.make_dup_gene_table(gene_names=DUP_GENES_FOR_TABLE,
                                regions=input.regions,
                                pb_SNV_stats=input.pb_SNV_stats,
                                il_SNV_stats=input.il_SNV_stats,
                                shared_SNV_stats=input.shared_SNV_stats,
                                pb_frac_map=input.pb_frac_map,
                                il_frac_map=input.il_frac_map,
                                output_file=output.tex)

rule filter_SNVs:
    params: job_name = 'filter_SNVs.{id}.{build}.{gene_name}.pb{pcov}x'
    input:  region = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/region.txt',
            pb_SNVs = 'data/{id}.{build}/variants/longshot.pacbio.blasr.{pcov}x.-z/all.GQ{pcov}.PASS.SNPs_ONLY.vcf.gz',
            il_SNVs = 'data/{id}.{build}/variants/illumina_30x.filtered/all.GQ50.PASS.SNPs_ONLY.vcf.gz',
            genome_file = 'genome_tracks/{build}.chrom.sizes.natural_order.txt'
    output: pb_SNVs = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/pacbio_SNVs.vcf.gz',
            il_SNVs = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/illumina_SNVs.vcf.gz',
            shared_SNVs = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/illumina_intersect_pacbio_SNVs.vcf.gz'
    shell:
        '''
        {RTGTOOLS} vcffilter --region=$(cat {input.region}) -i {input.pb_SNVs} -o {output.pb_SNVs}
        {RTGTOOLS} vcffilter --region=$(cat {input.region}) -i {input.il_SNVs} -o {output.il_SNVs}
        {BEDTOOLS} intersect -g {input.genome_file} -sorted -wa -header -a {output.pb_SNVs} -b {output.il_SNVs} | bgzip -c > {output.shared_SNVs}
        '''

rule gene_mappability:
    params: job_name = 'gene_mappability.{id}.{build}.{gene_name}.pb{pcov}x'
    input: region = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/region.txt',
           pb_bam = 'data/{id}.{build}/aligned_reads/pacbio/pacbio.blasr.all.{pcov}x.bam',
           pb_bai = 'data/{id}.{build}/aligned_reads/pacbio/pacbio.blasr.all.{pcov}x.bam.bai',
           il_bam = 'data/{id}.{build}/aligned_reads/illumina/illumina.30x.bam',
           il_bai = 'data/{id}.{build}/aligned_reads/illumina/illumina.30x.bam.bai',
           hg19    = 'data/genomes/hg19.fa',
           hg19_ix = 'data/genomes/hg19.fa.fai',
           hs37d5    = 'data/genomes/hs37d5.fa',
           hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
           hg38 = 'data/genomes/hg38.fa',
           hg38_ix = 'data/genomes/hg38.fa.fai'
    output: pb_num_map = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/pacbio_num_mappable.txt',
            pb_frac_map = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/pacbio_frac_mappable.txt',
            il_num_map = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/illumina_num_mappable.txt',
            il_frac_map = 'data/{id}.{build}/gene_mappability_pb{pcov}x/{gene_name}/illumina_frac_mappable.txt',
    run:
        (chrom, start, end) = parse_region_file(input.region)
        chrom_prefix = ''

        if wildcards.id == 'NA12878':
            chrom_prefix = 'chr'
            pb_ref_fa = input.hg19
            il_ref_fa = input.hs37d5
        elif wildcards.id == 'simulation':
            pb_ref_fa = input.hs37d5
            il_ref_fa = input.hs37d5
        else:
            pb_ref_fa = input.hg38
            il_ref_fa = input.hg38

        shell('''
        {MAP_COUNTER} -m --chrom {chrom_prefix}{chrom}:{start}-{end} --bam {input.pb_bam} --ref {pb_ref_fa} --min_cov 10 \
        --min_mapq 30 > {output.pb_num_map}
        ''')
        with open(output.pb_num_map, 'r') as inf, open(output.pb_frac_map, 'w') as outf:
            pb_num_map = float(inf.read().strip())
            pb_frac_map = pb_num_map / (end - start)
            print(pb_frac_map,end='',file=outf)

        shell('''
        {MAP_COUNTER} -m --chrom {chrom}:{start}-{end} --bam {input.il_bam} --ref {il_ref_fa} --min_cov 10 \
        --min_mapq 30 > {output.il_num_map}
        ''')
        with open(output.il_num_map, 'r') as inf, open(output.il_frac_map, 'w') as outf:
            il_num_map = float(inf.read().strip())
            il_frac_map = il_num_map / (end - start)
            print(il_frac_map,end='',file=outf)

rule get_region:
    params: job_name = 'get_region.{id}.{build}.pb{pcov}x'
    input:  gencode_bed = 'genome_tracks/gencode_{build}.bed'
    output: region = 'data/{id}.{build,(hg19|hg38)}/gene_mappability_pb{pcov}x/{gene_name}/region.txt'
    run:
        (chrom, start, end) = search_gencode_bed(wildcards.gene_name, input.gencode_bed)
        region_str = "{}:{}-{}".format(chrom, start+1, end) # convert to 1-based, inclusive string

        with open(output.region,'w') as outf:
            print(region_str,end="",file=outf)

rule get_region_1000g:
    params: job_name = 'get_region_1000g.{id}.pb{pcov}x'
    input:  gencode_bed = 'genome_tracks/gencode_hg19.bed'
    output: region = 'data/{id}.1000g/gene_mappability_pb{pcov}x/{gene_name}/region.txt'
    run:
        (chrom, start, end) = search_gencode_bed(wildcards.gene_name, input.gencode_bed)
        region_str = "{}:{}-{}".format(chrom, start+1, end) # convert to 1-based, inclusive string
        assert(region_str[:3] == 'chr')
        region_str = region_str[3:]
        with open(output.region,'w') as outf:
            print(region_str,end="",file=outf)

rule download_gencode_bed:
    params: job_name = 'download_gencode_bed.{build}'
    output: bed = 'genome_tracks/gencode_{build,(hg19|hg38)}.bed'
    run:
        release = 28 if wildcards.build == 'hg38' else 19
        gencode_url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{0}/gencode.v{0}.annotation.gff3.gz'.format(release)
        shell('''
        wget -qO- {gencode_url} \
            | gunzip --stdout - \
            | awk '$3 == "gene"' - \
            | convert2bed -i gff - \
            > {output.bed}
        ''') # taken from https://bioinformatics.stackexchange.com/questions/895/how-to-obtain-bed-file-with-coordinates-of-all-genes
