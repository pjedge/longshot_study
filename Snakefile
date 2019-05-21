import paper_tables_and_figures as ptf
from analyze_variants import analyze_variants, generate_random_calls
import time
from replace_empty_gt_with_reference import replace_empty_gt_with_reference
import sys
sys.path.append('HapCUT2/utilities')
import calculate_haplotype_statistics as chs
import prune_haplotype as ph
import pickle
import datetime
import pysam
from calculate_coverage import calculate_coverage
from math import sqrt

# DATA URLs
HG19_URL     = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit'
HS37D5_URL     = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
HG38_URL = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz' #'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
HG19_SDF_URL   = 'https://s3.amazonaws.com/rtg-datasets/references/hg19.sdf.zip'
TG_v37_SDF_URL = 'https://s3.amazonaws.com/rtg-datasets/references/1000g_v37_phase2.sdf.zip'
GRCh38_SDF_URL = 'https://s3.amazonaws.com/rtg-datasets/references/GRCh38.sdf.zip'
KNOWN_INDELS_URL_1000g = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz'
KNOWN_INDELS_URL_HG38 = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'

# PATHS TO TOOLS
TWOBITTOFASTA = 'twoBitToFa'
SAMTOOLS       = 'samtools'
FASTQ_DUMP     = 'fastq-dump'
LONGSHOT_BETA_ALN_PARAMS  = 'longshot_beta_const_aln_params/target/release/longshot'
LONGSHOT_BETA  = 'longshot_beta/target/release/longshot'
LONGSHOT       = 'longshot'
RTGTOOLS       = 'rtg'
BGZIP = 'bgzip'
TABIX = 'tabix'
FREEBAYES      = '/home/pedge/git/freebayes/bin/freebayes' # v1.0.2-33-gdbb6160
WHATSHAP = 'whatshap'
SIMLORD = 'simlord'
DWGSIM = 'dwgsim'
BWA = 'bwa'
BLASR = 'blasr'
BAX2BAM = 'bax2bam'
BAM2FASTQ = 'bam2fastq'
SAWRITER = 'sawriter'
NGMLR = 'ngmlr'
MINIMAP2 = 'minimap2'
BCFTOOLS = 'bcftools'
PYFAIDX = 'faidx'
BEDTOOLS = 'bedtools' # v 2.27
EXTRACTHAIRS = 'HapCUT2/build/extractHAIRS'
HAPCUT2 = 'HapCUT2/build/HAPCUT2'
PYTHON = 'python3'
MAP_COUNTER = 'map_counter/target/release/map_counter'
SEQTK = 'seqtk'

chroms = ['{}'.format(i) for i in range(1,23)]
ref_file = {'1000g':'data/genomes/hs37d5.fa', 'hg38':'data/genomes/hg38.fa'}

include: "simulation.1000g.snakefile"
#include: "simulation.test.snakefile"
include: "NA12878.1000g.snakefile"
include: "NA12878.hg38.snakefile"
include: "NA24385.hg38.snakefile"  # AJ Son,    hg38
include: "NA24385.1000g.snakefile"  # AJ Son,   1000g
include: "NA24143.hg38.snakefile"  # AJ Mother, hg38
include: "NA24143.1000g.snakefile"  # AJ Mother, 1000g
include: "NA24149.hg38.snakefile"  # AJ Father, hg38
include: "NA24149.1000g.snakefile"  # AJ Father, 1000g
include: "aj_trio.snakefile" #
include: "paper_tables_and_figures.snakefile"
include: "haplotyping.snakefile"
include: "make_variant_counts_table.snakefile"
include: "duplicated_gene_visualization.snakefile"
include: "make_duplicated_gene_table.snakefile"
include: "map_giab_reads.snakefile"

#rule default:
 #  input:
#        'data/plots/longshot_params_comparison_30x.png'
        #'data/plots/haplotyping_results_barplot.png',
        #'data/plots/fig3_precision_recall_bars_NA12878_AJ_Trio_with_haplotyping_results.blasr.hg38.png',
        #'data/NA24385.1000g/vcfeval/clairvoyante.pacbio.ngmlr.69x.unfiltered/all',
        #'data/NA24385.1000g/vcfeval/whatshap.pacbio.ngmlr.69x.dp_dn_filtered/all',
        #'data/NA24385.1000g/vcfeval/longshot.pacbio.ngmlr.69x._/all',
        #'data/NA12878.1000g/vcfeval/clairvoyante.pacbio.ngmlr.44x.unfiltered/all',
        #'data/NA12878.1000g/vcfeval/whatshap.pacbio.ngmlr.44x.dp_dn_filtered/all',
        #'data/NA12878.1000g/vcfeval/longshot.pacbio.ngmlr.44x._/all',
        #'data/output/longshot_v_whatshap.aj_son_hg38_blasr.all.tex'
        #'data/NA12878.hg38/vcfeval/whatshap.ont.minimap2.30x.unfiltered/all',
        #'data/NA24385.hg38/vcfeval/whatshap.pacbio.blasr.69x.dp_dn_filtered/all',
        #expand('data/NA24385.hg38/vcfeval/clairvoyante.pacbio.blasr.69x.unfiltered/{c}',c=chr20to22),
        #expand('data/NA12878.1000g/vcfeval/clairvoyante.pacbio.blasr.44x.unfiltered/{c}',c=chr20to22),
        #expand('data/NA24385.hg38/variants/genomicconsensus.pacbio.blasr.69x.no_sample_info/{c}.vcf',c=chr20to22),
        #'data/NA24385.hg38/variants/genomicconsensus.pacbio.blasr.69x.no_sample_info/all.vcf',
        #'data/NA12878.1000g/vcfeval/genomicconsensus.pacbio.blasr.44x.unfiltered/all'
        #'data/NA12878.1000g/vcfeval/clairvoyante.pacbio.blasr.44x.unfiltered/20',
        #'data/NA12878.1000g/vcfeval/clairvoyante.pacbio.blasr.44x.unfiltered/5',
        #'data/NA12878.1000g/vcfeval/clairvoyante.pacbio.blasr.44x.unfiltered/4',

        #'data/NA12878.1000g/vcfeval/clairvoyante.pacbio.blasr.44x.unfiltered/3',

        #'data/NA12878.1000g/vcfeval/clairvoyante.pacbio.blasr.44x.unfiltered/all',
        #'data/NA12878.1000g/vcfeval/whatshap.pacbio.blasr.44x.unfiltered/all',
        #'data/NA12878.1000g/vcfeval/whatshap.pacbio.blasr.44x.dp_filtered/all',
        #'data/NA12878.1000g/vcfeval/whatshap.pacbio.blasr.44x.dp_dn_filtered/all',

#        'data/NA12878.hg38/vcfeval/whatshap.ont.minimap2.30x._/all'


rule all:
    input:
        # tables & figures
        'data/plots/3method_NA12878.1000g.ngmlr.prec_recall_all.png',
        'data/plots/longshot_params_comparison_44x.png',
        'data/plots/NA12878.hg38_ONT_PR_curve_all.png',
        'data/plots/fig3_precision_recall_bars_NA12878_AJ_Trio_with_haplotyping_results.blasr.hg38.png',
        'data/plots/prec_recall_4panel_blasr.all.png',
        'data/plots/actual_vs_effective_coverage.chr1.NA12878.44x.png',
        'data/output/prec_recall_table_known_indels_filtered.tex',
        'data/plots/NA12878_variants_outside_GIAB_confident_venn_diagram.png',
        'data/output/four_GIAB_genomes_table_extended.aj_trio_hg38_blasr.all.tex',     # table 1
        'data/output/variant_analysis_fp_fn_NA12878.1000g.blasr.44x.GQ44.1.tex',       # table 2
        'data/output/variant_counts_table.NA12878.1000g.il30x.blasr.pb30x.GQ30.tex',   # table 3
        'data/output/pacbio_mendelian_concordance_table.blasr.tex',           # table 4
        'data/plots/haplotyping_results_barplot.png',
        'data/plots/plot_mappability_bars.simulation.1000g.png',
        'data/plots/simulation_pr_barplot_genome_vs_segdup.all.GQ50.png',
        'data/plots/simulation_pr_barplot_genome_vs_segdup_extended.all.GQ50.png',
        'data/plots/NA12878_variants_outside_GIAB_confident_inside_PG_confident_venn_diagram.png',\
        'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.read_lengths',
        'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.read_lengths',
        'data/output/NA12878_chr1_45x_haplotype_assigned_N50_analysis.txt',
        #'data/NA12878.1000g/variants/misc/INTERSECT_PG_longshot.pacbio.blasr.44x._.all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.stats',
        #'data/NA12878.1000g/variants/misc/MINUS_longshot.pacbio.blasr.44x._.PG.all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.stats',
        #'data/NA12878.1000g/variants/misc/MINUS_PG_longshot.pacbio.blasr.44x._.all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.stats',


rule coding_exons_segdup_NA12878:
    params: job_name = 'coding_exons_segdup_NA12878.44x',
    input:  bam = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam',
            bai = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.bai',
            ref = 'data/genomes/hg19.fa',
            bed = 'genome_tracks/coding_exons_intersect_segmental_duplications_1000g_0.95_similar.bed'
    output: txt = 'data/output/NA12878.1000g_num_bases_coding_exons_covered.txt',
    run:
        shell('''echo '' > {output}''')
        with open(input.bed, 'r') as inf:
            for line in inf:
                el = line.strip().split()
                reg = 'chr'+el[0]+':'+str(int(el[1])+1)+'-'+el[2] # convert 0-index exclusive to 1-index inclusive
                shell('''
                {MAP_COUNTER} --chrom {reg} --bam {input.bam} --ref {input.ref} --min_cov 20 \
                --min_mapq 30 --map_frac 0.9 >> {output}
                ''')

rule run_pileups:
    params: job_name = 'run_pileups_{individual}.{cov}x',
    input:  bam = 'data/{individual}.hg38/aligned_reads/pacbio/pacbio.blasr.all.{cov}x.bam',
            pos = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/all.positions_only.txt'
    output: txt = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/{individual}.{cov}x.pileup_over_mendelian_positions.txt',
    shell: 'python3 run_pileups.py {input.bam} {input.pos} > {output.txt}'

# this function reads in a file containing a single integer, and returns that integer.
# this is designed for bam.median_coverage files that accompany a bam file and
# tell us the median coverage we've measured for that bam.
def parse_int_file(filename):
    with open(filename,'r') as inf:
        return int(float(inf.read().strip()))

# this function takes in a chromosome name in 1..22,X and a "genome build" name
# and appends a 'chr' prefix to the chrom name if the genome build is hg38
def chr_prefix(chrom, build):
    assert(build in ['1000g','hg38'])
    assert(chrom in chroms)
    if build == 'hg38':
        return 'chr'+chrom
    else:
        return chrom

def seconds_to_formatted_time(ss_time):
    hh = int(ss_time / 3600)
    ss_time -= float(hh*3600)
    mm = int(ss_time / 60)
    ss_time -= float(mm*60)
    ss = int(ss_time)
    return "{}:{:02d}:{:02d}".format(hh,mm,ss)

def formatted_time_to_seconds(fmt_time):
    vals = [float(val) for val in fmt_time.strip().split(':')]
    assert(len(vals) == 3)
    hh, mm, ss = vals
    return hh*3600 + mm*60 + ss

rule get_read_lengths:
    params: job_name = lambda wildcards: 'get_read_lengths.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.bam'
    output: '{x}.bam.read_lengths'
    shell: '''{SAMTOOLS} view {input} | awk '{{print length($10)}}' > {output}'''

rule remove_x_y_region_filter:
    params: job_name = 'remove_x_y_region_filter'
    input:  'data/{individual}.{build}/variants/ground_truth/region_filter.bed'
    output: 'data/{individual}.{build}/variants/ground_truth/region_filter.1_22_only.bed'
    shell:  'cat {input} | python3 filter_bed_chroms.py > {output}'

# rule vcfeval_rtgtools:
#     params: job_name = 'vcfeval_rtgtools.{individual}.{build}.{calls_name}.{chrom}',
#             region_arg = lambda wildcards: '--region={}'.format(chr_prefix(wildcards.chrom,wildcards.build)) if wildcards.chrom != 'all' else ''
#     input:  calls_vcf = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf.gz',
#             calls_ix = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf.gz.tbi',
#             ground_truth = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz',
#             ground_truth_ix = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz.tbi',
#             region_filter ='data/{individual}.{build}/variants/ground_truth/region_filter.1_22_only.bed',
#             sdf = 'data/genomes/{build}.sdf'
#     output: dir = directory('data/{individual}.{build}/vcfeval/{calls_name}/{chrom}')
#     shell:
#         '''
#         {RTGTOOLS} RTG_MEM=12g vcfeval \
#         {params.region_arg} \
#         -c {input.calls_vcf} \
#         -b {input.ground_truth} \
#         -e {input.region_filter} \
#         -t {input.sdf} \
#         -o {output.dir}
#         '''

rule vcfeval_rtgtools:
    params: job_name = 'vcfeval_rtgtools.{individual}.{build}.{calls_name}.{chrom}',
            region_arg = lambda wildcards: '--region={}'.format(chr_prefix(wildcards.chrom,wildcards.build)) if wildcards.chrom != 'all' else ''
    input:  calls_vcf = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf.gz',
            calls_ix = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf.gz.tbi',
            ground_truth = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz',
            ground_truth_ix = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz.tbi',
            region_filter ='data/{individual}.{build}/variants/ground_truth/region_filter.1_22_only.bed',
            sdf = 'data/genomes/{build}.sdf'
    output: dir = directory('data/{individual}.{build}/vcfeval/{calls_name}/{chrom}')
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcfeval \
        {params.region_arg} \
        -c {input.calls_vcf} \
        -b {input.ground_truth} \
        -e {input.region_filter} \
        -t {input.sdf} \
        --decompose \
        -o {output.dir}
        '''

rule vcfeval_rtgtools_last3:
    params: job_name = 'vcfeval_rtgtools.{individual}.{build}.{calls_name}.all',
    input:  calls_vcf = 'data/{individual}.{build}/variants/{calls_name}/all.vcf.gz',
            calls_ix = 'data/{individual}.{build}/variants/{calls_name}/all.vcf.gz.tbi',
            ground_truth = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz',
            ground_truth_ix = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz.tbi',
            region_filter ='data/{individual}.{build}/variants/ground_truth/region_filter.last3.bed',
            sdf = 'data/genomes/{build}.sdf'
    output: dir = directory('data/{individual}.{build}/vcfeval/{calls_name}/last3')
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcfeval \
        -c {input.calls_vcf} \
        -b {input.ground_truth} \
        -e {input.region_filter} \
        -t {input.sdf} \
        --decompose \
        -o {output.dir}
        '''

rule rtg_filter_SNVs_known_indels:
    params: job_name = 'rtg_filter_SNVs_known_indels.{individual}.{build}.longshot.pacbio.blasr.{cov}x',
    input:  vcfgz = 'data/{individual}.{build}/variants/longshot.pacbio.blasr.{cov}x._/all.vcf.gz',
            tbi = 'data/{individual}.{build}/variants/longshot.pacbio.blasr.{cov}x._/all.vcf.gz.tbi',
            known_indels = 'genome_tracks/Mills_and_1000G_gold_standard.indels.{build}.10bp_window.bed'
    output: vcf = 'data/{individual}.{build}/variants/known_indels_removed.longshot.pacbio.blasr.{cov}x._/all.vcf',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-bed {input.known_indels} -i {input.vcfgz} -o {output.vcf}.gz; gunzip {output.vcf}.gz'

from intersect_SNV_genotypes import intersect_SNV_genotypes
rule intersect_longshot_clairvoyante:
    params: job_name = 'intersect_longshot_clairvoyante.{individual}.{build}.{cov}',
    input:  ls  = 'data/{individual}.{build}/variants/longshot.pacbio.ngmlr.{cov}x._/all.vcf',
            clv = 'data/{individual}.{build}/variants/clairvoyante.pacbio.ngmlr.{cov}x.unfiltered/all.vcf'
    output: vcf = 'data/{individual}.{build}/variants/longshot_x_clairvoyante.longshot.pacbio.ngmlr.{cov}x._/all.vcf',
    run:
        intersect_SNV_genotypes(input.clv, input.ls, output.vcf)

from intersect_SNV_genotypes import intersect_SNV_genotypes
rule intersect_longshot_clairvoyante_blasr_ngmlr:
    params: job_name = 'intersect_longshot_clairvoyante_blasr_ngmlr.{individual}.{build}.{cov}',
    input:  ls  = 'data/{individual}.{build}/variants/longshot.pacbio.blasr.{cov}x._/all.vcf',
            clv = 'data/{individual}.{build}/variants/clairvoyante.pacbio.ngmlr.{cov}x.unfiltered/all.vcf'
    output: vcf = 'data/{individual}.{build}/variants/longshot_x_clairvoyante.longshot.pacbio.blasr_ngmlr.{cov}x._/all.vcf',
    run:
        intersect_SNV_genotypes(input.clv, input.ls, output.vcf)

def convert_known_indels_to_bed(input, output, bp_pad=5):
    with pysam.VariantFile(input) as vcf, open(output,'w') as outf:
        for record in vcf:
            alt_len = max([len(x) for x in record.alleles])
            l_pos = record.pos-bp_pad
            r_pos = record.pos + alt_len + bp_pad
            print("{}\t{}\t{}".format(record.chrom, l_pos, r_pos),file=outf)

rule convert_known_indels_to_bed:
    params: job_name = 'convert_known_indels_{build}_to_bed',
    input:  'genome_tracks/Mills_and_1000G_gold_standard.indels.{build}.vcf.gz'
    output: 'genome_tracks/Mills_and_1000G_gold_standard.indels.{build}.10bp_window.bed'
    run:
        convert_known_indels_to_bed(input[0], output[0], bp_pad=5)

rule download_known_indels_1000g:
    params: job_name = 'download_known_indels_1000g',
    output: vcfgz_bad = 'genome_tracks/Mills_and_1000G_gold_standard.indels.1000g.needs_recompression.vcf.gz',
            vcfgz = 'genome_tracks/Mills_and_1000G_gold_standard.indels.1000g.vcf.gz'
    shell:
        '''
        wget {KNOWN_INDELS_URL_1000g} -O {output.vcfgz_bad}
        gunzip -c {output.vcfgz_bad} | bgzip -c > {output.vcfgz}
        '''

rule download_known_indels_hg38:
    params: job_name = 'download_known_indels_hg38',
    output: vcfgz = 'genome_tracks/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
    shell: 'wget {KNOWN_INDELS_URL_HG38} -O {output.vcfgz}'

rule rtg_filter_SNVs_ground_truth:
    params: job_name = 'rtg_filter_SNVs_ground_truth.{individual}.{build}',
    input:  vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.vcf.gz.tbi'
    output: vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --snps-only -i {input.vcfgz} -o {output.vcfgz}'

rule rtg_decompose_variants_ground_truth:
    params: job_name = 'rtg_decompose_variants_ground_truth.{individual}.{build}',
    input:  vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz.tbi'
    output: vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfdecompose --break-mnps --break-indels -i {input.vcfgz} -o {output.vcfgz}'

rule filter_illumina_SNVs:
    params: job_name = 'filter_SNVs_illumina.{individual}.{build}.chr{chrom}',
    input:  vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov}x.unfiltered/{chrom}.vcf.gz',
            runtime = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov}x.unfiltered/{chrom}.vcf.runtime',
            cov = 'data/{individual}.{build}/aligned_reads/illumina/illumina.aligned.all.{cov}x.bam.median_coverage'
    output: vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov,\d+}x.filtered/{chrom,(\d+)}.vcf',
            runtime = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov,\d+}x.filtered/{chrom,(\d+)}.vcf.runtime'
    run:
        median_cov = parse_int_file(input.cov)
        max_cov = int(median_cov + 5*sqrt(median_cov))
        shell('{RTGTOOLS} RTG_MEM=12g vcffilter -D {max_cov} -i {input.vcfgz} -o {output.vcf}.gz')
        shell('gunzip -c {output.vcf}.gz > {output.vcf}')
        shell('cp {input.runtime} {output.runtime}')

rule combine_chrom:
    params: job_name = 'combine_chroms.{individual}.{build}.{caller}.{tech}.{further_info}',
    input: expand('data/{{individual}}.{{build}}/variants/{{caller}}.{{tech}}.{{further_info}}/{chrom}.vcf',chrom=chroms)
    output: 'data/{individual}.{build}/variants/{caller,(freebayes|longshot|whatshap|genomicconsensus|clairvoyante)}.{tech,(pacbio|illumina|ont)}.{further_info}/all.vcf'
    shell:
        '''
        grep -P '^#' {input[0]} > {output}; # grep header
        cat {input} | grep -Pv '^#' >> {output}; # cat files, removing the headers.
        '''

hg19_chroms = set(['chr{}'.format(i) for i in range(1,23)])
hs37d5_chroms = set([str(i) for i in range(1,23)])
def remove_chr_from_vcf(in_vcf, out_vcf):
    with open(in_vcf, 'r') as inf, open(out_vcf, 'w') as outf:
        for line in inf:
            if line[0] == '#':
                print(line.strip(),file=outf)
                continue
            el = line.strip().split('\t')
            assert(el[0] in hg19_chroms)
            el[0] = el[0][3:]
            assert(el[0] in hs37d5_chroms)
            print("\t".join(el),file=outf)

rule add_runtimes:
    params: job_name = 'add_runtimes.{individual}.{build}.{calls_name}',
    input: expand('data/{{individual}}.{{build}}/variants/{{calls_name}}/{chrom}.vcf.runtime',chrom=chroms)
    output: 'data/{individual}.{build}/variants/{calls_name}/all.vcf.runtime'
    run:
        t = 0 # time in seconds
        for f in input:
            with open(f,'r') as inf:
                t += formatted_time_to_seconds(inf.readline().strip())

        with open(output[0],'w') as outf:
            print(seconds_to_formatted_time(t),file=outf)

rule run_longshot_pacbio:
    params: job_name = 'longshot.pacbio.{aligner}.{individual}.{build}.cov{cov}.{options}.chr{chrom}',
            debug = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.debug'
    input:  bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam.bai',
            cov = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam.median_coverage',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.vcf',
            potential_snvs = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.debug/1.0.potential_SNVs.vcf',
            runtime = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{cov}x.{options}/{chrom,(\d+|X|Y)}.vcf.runtime'
    run:
        median_cov = parse_int_file(input.cov)
        max_cov = int(median_cov + 5*sqrt(median_cov))
        options_str = wildcards.options.replace('_',' ')
        if wildcards.individual == 'NA12878' and wildcards.build == '1000g' and wildcards.aligner == 'blasr':
            t1 = time.time()
            shell('{LONGSHOT} -r chr{wildcards.chrom} -F -C {max_cov} -d {params.debug} {options_str} -s {wildcards.individual} --bam {input.bam} --ref {input.hg19} --out {output.vcf}.tmp')
            t2 = time.time()
            # remove 'chr' from reference name in vcf
            remove_chr_from_vcf(output.vcf+'.tmp',output.vcf)
            # remove 'chr' from no-haplotype version of vcf
        else:
            w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
            w_ref = ref_file[wildcards.build]
            t1 = time.time()
            shell('{LONGSHOT} -r {w_chrom} -F -C {max_cov} -d {params.debug} {options_str} -s {wildcards.individual} --bam {input.bam} --ref {w_ref} --out {output.vcf}')
            t2 = time.time()

        with open(output.runtime,'w') as outf:
            print(seconds_to_formatted_time(t2-t1),file=outf)

rule run_longshot_beta_pacbio:
    params: job_name = 'longshot_beta.pacbio.{aligner}.{individual}.{build}.cov{cov}.{options}.chr{chrom}',
    input:  bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam.bai',
            cov = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam.median_coverage',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/longshot_beta.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.vcf',
            debug = directory('data/{individual}.{build}/variants/longshot_beta.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.debug'),
            runtime = 'data/{individual}.{build}/variants/longshot_beta.pacbio.{aligner}.{cov}x.{options}/{chrom,(\d+|X|Y)}.vcf.runtime'
    run:
        median_cov = parse_int_file(input.cov)
        max_cov = int(median_cov + 5*sqrt(median_cov))
        options_str = wildcards.options.replace('_',' ')
        if wildcards.individual == 'NA12878' and wildcards.build == '1000g' and wildcards.aligner == 'blasr':
            t1 = time.time()
            shell('{LONGSHOT_BETA} -r chr{wildcards.chrom} -F -C {max_cov} -d {output.debug} {options_str} -s {wildcards.individual} --bam {input.bam} --ref {input.hg19} --out {output.vcf}.tmp')
            t2 = time.time()
            # remove 'chr' from reference name in vcf
            remove_chr_from_vcf(output.vcf+'.tmp',output.vcf)
            # remove 'chr' from no-haplotype version of vcf
        else:
            w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
            w_ref = ref_file[wildcards.build]
            t1 = time.time()
            shell('{LONGSHOT_BETA} -r {w_chrom} -F -C {max_cov} -d {output.debug} {options_str} -s {wildcards.individual} --bam {input.bam} --ref {w_ref} --out {output.vcf}')
            t2 = time.time()

        with open(output.runtime,'w') as outf:
            print(seconds_to_formatted_time(t2-t1),file=outf)

rule run_longshot_ont:
    params: job_name = 'longshot.ont.{aligner}.{individual}.{build}.cov{cov}.{options}.chr{chrom}',
            debug = 'data/{individual}.{build}/variants/longshot.ont.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.debug'
    input:  bam = 'data/{individual}.{build}/aligned_reads/ont/ont.{aligner}.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/ont/ont.{aligner}.all.{cov}x.bam.bai',
            cov = 'data/{individual}.{build}/aligned_reads/ont/ont.{aligner}.all.{cov}x.bam.median_coverage',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/longshot.ont.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.vcf',
            potential_snvs = 'data/{individual}.{build}/variants/longshot.ont.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+|X|Y)}.debug/1.0.potential_SNVs.vcf',
            runtime = 'data/{individual}.{build}/variants/longshot.ont.{aligner}.{cov}x.{options}/{chrom,(\d+|X|Y)}.vcf.runtime'
    run:
        median_cov = parse_int_file(input.cov)
        max_cov = int(median_cov + 5*sqrt(median_cov))
        options_str = wildcards.options.replace('_',' ')
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        w_ref = ref_file[wildcards.build]
        t1 = time.time()
        shell('{LONGSHOT} -r {w_chrom} -F -C {max_cov} -d {params.debug} {options_str} -s {wildcards.individual} --bam {input.bam} --ref {w_ref} --out {output.vcf}')
        t2 = time.time()

        with open(output.runtime,'w') as outf:
            print(seconds_to_formatted_time(t2-t1),file=outf)

from filter_SNVs import filter_SNVs_density
rule filter_othermethod_SNVs_depth:
    params: job_name = 'filter_SNVs_{method}_depth.{individual}.{build}.chr{chrom}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{method}.{tech}.{aligner}.{cov}x.unfiltered/{chrom}.vcf.gz',
            runtime = 'data/{individual}.{build}/variants/{method}.{tech}.{aligner}.{cov}x.unfiltered/{chrom}.vcf.runtime',
            cov = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam.median_coverage'
    output: vcf = 'data/{individual}.{build}/variants/{method}.{tech,(pacbio|ont)}.{aligner}.{cov,\d+}x.dp_filtered/{chrom,(\d+)}.vcf',
            runtime = 'data/{individual}.{build}/variants/{method}.{tech}.{aligner}.{cov,\d+}x.dp_filtered/{chrom,(\d+)}.vcf.runtime'
    run:
        median_cov = parse_int_file(input.cov)
        max_cov = int(median_cov + 5*sqrt(median_cov))
        shell('{RTGTOOLS} RTG_MEM=12g vcffilter -D {max_cov} -i {input.vcfgz} -o {output.vcf}.gz')
        shell('gunzip -c {output.vcf}.gz > {output.vcf}')
        shell('cp {input.runtime} {output.runtime}')

rule filter_othermethod_SNVs_density:
    params: job_name = 'filter_SNVs_{method}_density.{individual}.{build}.chr{chrom}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{method}.{tech}.{aligner}.{cov}x.dp_filtered/{chrom}.vcf',
            runtime = 'data/{individual}.{build}/variants/{method}.{tech}.{aligner}.{cov}x.dp_filtered/{chrom}.vcf.runtime',
    output: vcf = 'data/{individual}.{build}/variants/{method}.{tech,(pacbio|ont)}.{aligner}.{cov,\d+}x.dp_dn_filtered/{chrom,(\d+)}.vcf',
            runtime = 'data/{individual}.{build}/variants/{method}.{tech}.{aligner}.{cov,\d+}x.dp_dn_filtered/{chrom,(\d+)}.vcf.runtime'
    run:
        filter_SNVs_density(input.vcfgz, output.vcf, density_count=10, density_len=500, density_qual=50)
        shell('cp {input.runtime} {output.runtime}')

rule run_whatshap_pacbio:
    params: job_name = 'whatshap.pacbio.{aligner}.{individual}.{build}.cov{cov}.chr{chrom}',
    input:  bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{cov}x.bam.bai',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai',
            potential_snvs = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{cov}x._/{chrom}.debug/1.0.potential_SNVs.vcf',
    output: vcf = 'data/{individual}.{build}/variants/whatshap.pacbio.{aligner}.{cov,\d+}x.unfiltered/{chrom,(\d+|X|Y)}.vcf',
            runtime = 'data/{individual}.{build}/variants/whatshap.pacbio.{aligner}.{cov}x.unfiltered/{chrom,(\d+|X|Y)}.vcf.runtime'
    run:
        if wildcards.individual == 'NA12878' and wildcards.build == '1000g' and wildcards.aligner == 'blasr':
            t1 = time.time()
            shell('{WHATSHAP} genotype --ignore-read-groups --reference {input.hg19} -o {output.vcf}.tmp {input.potential_snvs} {input.bam}')
            t2 = time.time()
            # remove 'chr' from reference name in vcf
            remove_chr_from_vcf(output.vcf+'.tmp',output.vcf)
            # remove 'chr' from no-haplotype version of vcf
        else:
            w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
            w_ref = ref_file[wildcards.build]
            t1 = time.time()
            shell('{WHATSHAP} genotype --ignore-read-groups --reference {w_ref} -o {output.vcf} {input.potential_snvs} {input.bam}')
            t2 = time.time()

        with open(output.runtime,'w') as outf:
            print(seconds_to_formatted_time(t2-t1),file=outf)

rule run_whatshap_ont:
    params: job_name = 'whatshap.ont.{aligner}.{individual}.{build}.cov{cov}.chr{chrom}',
    input:  bam = 'data/{individual}.{build}/aligned_reads/ont/ont.{aligner}.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/ont/ont.{aligner}.all.{cov}x.bam.bai',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai',
            potential_snvs = 'data/{individual}.{build}/variants/longshot.ont.{aligner}.{cov}x._/{chrom}.debug/1.0.potential_SNVs.vcf',
    output: vcf = 'data/{individual}.{build}/variants/whatshap.ont.{aligner}.{cov,\d+}x.unfiltered/{chrom,(\d+|X|Y)}.vcf',
            runtime = 'data/{individual}.{build}/variants/whatshap.ont.{aligner}.{cov}x.unfiltered/{chrom,(\d+|X|Y)}.vcf.runtime'
    run:
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        w_ref = ref_file[wildcards.build]
        t1 = time.time()
        shell('{WHATSHAP} genotype --ignore-read-groups --reference {w_ref} -o {output.vcf} {input.potential_snvs} {input.bam}')
        t2 = time.time()

        with open(output.runtime,'w') as outf:
            print(seconds_to_formatted_time(t2-t1),file=outf)

# before this rule is executed, it is necessary to modify the conda env
# run snakemake with the command:
# snakemake --use-conda --create-envs-only
# this will create clairvoyante's conda environment
# then, find and activate that environment (should have no name and long hash directory path) e.g.:
# conda activate .snakemake/conda/7c00edba
# and run:
# wget https://bootstrap.pypa.io/get-pip.py
# pypy get-pip.py
# pypy -m pip install --no-cache-dir intervaltree==2.1.0
rule run_clairvoyante:
    params: job_name = 'clairvoyante.{tech}.{aligner}.{individual}.{build}.cov{cov}.chr{chrom}',
            w_chrom = lambda wildcards: chr_prefix(wildcards.chrom, wildcards.build),
            w_ref = lambda wildcards: ref_file[wildcards.build]
    input:  bam = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam.bai',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/clairvoyante.{tech,(pacbio|ont)}.{aligner}.{cov,\d+}x.unfiltered/{chrom,(\d+|X|Y)}.vcf',
            runtime = 'data/{individual}.{build}/variants/clairvoyante.{tech}.{aligner}.{cov}x.unfiltered/{chrom,(\d+|X|Y)}.vcf.runtime'
    conda:  'envs/clairvoyante.yaml'
    script: 'scripts/clairvoyante.py'

CLAIRVOYANTE_MODELS_URL = 'http://www.bio8.cs.hku.hk/trainedModels.tbz'
#
rule download_clairvoyante_models:
    params: job_name = 'download_clairvoyante_models',
    output: dir('data/clairvoyante_models')
    shell: 'curl {CLAIRVOYANTE_MODELS_URL} | tar -jxf -'

rule run_genomicconsensus:
    params: job_name = 'genomicconsensus.{tech}.{aligner}.{individual}.{build}.cov{cov}.chr{chrom}',
            w_chrom = lambda wildcards: chr_prefix(wildcards.chrom, wildcards.build),
            w_ref = lambda wildcards: ref_file[wildcards.build]
    input:  bam = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam.bai',
            pbi = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam.pbi',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/genomicconsensus.{tech,(pacbio|ont)}.{aligner}.{cov,\d+}x.no_sample_info/{chrom,(\d+|X|Y)}.vcf',
            runtime = 'data/{individual}.{build}/variants/genomicconsensus.{tech}.{aligner}.{cov}x.unfiltered/{chrom,(\d+|X|Y)}.vcf.runtime'
    conda:  'envs/genomicconsensus.yaml'
    script: 'scripts/genomicconsensus.py'

rule pbindex_bam:
    params: job_name = lambda wildcards: 'pbindex_bam.{}'.format(str(wildcards.x).replace("/", "."))
    input:  bam = '{x}.bam'
    output: bai = '{x}.bam.pbi'
    shell:  'pbindex {input.bam}'


# Call 30x Illumina variants
#rule call_variants_Illumina:
#    params: job_name = 'call_illumina.{individual}.{build}.{cov}x.chr{chrom}',
#    input: bam = 'data/{individual}.{build}/aligned_reads/illumina/illumina.aligned.all.{cov}x.bam',
#            bai = 'data/{individual}.{build}/aligned_reads/illumina/illumina.aligned.all.{cov}x.bam.bai',
#            ref_1000g_fa = 'data/genomes/hs37d5.fa',
#            ref_1000g_fai = 'data/genomes/hs37d5.fa.fai',
#            ref_hg38_fa = 'data/genomes/hg38.fa',
#            ref_hg38_fai = 'data/genomes/hg38.fa.fai'
#    output: vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov}x.unfiltered/{chrom,(\d+)}.vcf',
#            runtime = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov}x.unfiltered/{chrom,(\d+)}.vcf.runtime'
#    conda: 'envs/freebayes.yaml'
#    script: 'scripts/call_variants_Illumina.py'

# Call 30x Illumina variants
rule call_variants_Illumina:
    params: job_name = 'call_illumina.{individual}.{build}.{cov}x.chr{chrom}',
    input: bam = 'data/{individual}.{build}/aligned_reads/illumina/illumina.aligned.all.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/illumina/illumina.aligned.all.{cov}x.bam.bai',
            ref_1000g_fa = 'data/genomes/hs37d5.fa',
            ref_1000g_fai = 'data/genomes/hs37d5.fa.fai',
            ref_hg38_fa = 'data/genomes/hg38.fa',
            ref_hg38_fai = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov}x.unfiltered/{chrom,(\d+)}.vcf',
            runtime = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{cov}x.unfiltered/{chrom,(\d+)}.vcf.runtime'
    run:
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        w_ref = ref_file[wildcards.build]
        t1 = time.time()
        shell('''
        {FREEBAYES} -f {w_ref} \
        --standard-filters \
        --region {w_chrom} \
         --genotype-qualities \
         {input.bam} \
          > {output.vcf}
        ''')
        t2 = time.time()
        with open(output.runtime,'w') as outf:
            print(seconds_to_formatted_time(t2-t1),file=outf)

rule generate_genomecov_bed:
    params: job_name = 'generate_genomecov_bed.{individual}.{build}.{tech}.{info}.MAPQ{mapq}'
    input:  expand('data/{{individual}}.{{build}}/aligned_reads/{{tech}}/genomecov_histograms_mapq{{mapq}}/{{tech}}.{{aligner}}.30x.{{chrom}}.txt', chrom=chroms)
    output: 'data/{individual}.{build}/aligned_reads/{tech}/genomecov_histograms_mapq{mapq}/{tech}.{aligner}.30x.all.txt'
    shell: 'cat {input} > {output}'

rule generate_coverage_histogram:
    params: job_name = 'generate_coverage_histogram.{tech}.{individual}.{build}.MAPQ{mapq}.{chrom}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.{chrom}.30x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.{chrom}.30x.bam.bai'
    output: 'data/{individual}.{build}/aligned_reads/{tech}/genomecov_histograms_mapq{mapq}/{tech}.{aligner}.30x.{chrom,(\d+)}.txt'
    run:
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        # NA12878 is an exception, actually has hg19 chrom names instead of 1000g chrom names
        if wildcards.individual == 'NA12878' and wildcards.build == '1000g' and wildcards.tech=='pacbio':
            w_chrom = 'chr' + w_chrom
        shell('''
        {SAMTOOLS} view -F 3844 -q {wildcards.mapq} {input.bam} -hb {w_chrom} | \
        {BEDTOOLS} genomecov -ibam - | \
        grep -P '^{w_chrom}\\t' > {output}
        ''')

# download hg19 reference, for the aligned pacbio reads
rule download_hg19:
    params: job_name = 'download_hg19',
            url      = HG19_URL
    output: 'data/genomes/hg19.fa'
    shell:
        '''
        wget {params.url} -O {output}.2bit
        {TWOBITTOFASTA} {output}.2bit {output}
        '''

rule download_hg38:
    params: job_name = 'download_hg38',
    output: 'data/genomes/hg38.fa'
    shell: 'wget {HG38_URL} -O {output}.gz; gunzip {output}.gz'

rule download_hs37d5:
    params: job_name = 'download_hs37d5',
    output: 'data/genomes/hs37d5.fa'
    shell: 'wget {HS37D5_URL} -O {output}.gz; gunzip {output}.gz'

# download hg19 reference, for the aligned pacbio reads
rule download_hg19_sdf:
    params: job_name = 'download_hg19_sdf',
    output: directory('data/genomes/hg19.sdf')
    shell:
        '''
        wget {HG19_SDF_URL} -O {output}.zip;
        unzip {output}.zip -d data/genomes
        '''

# download 1000g_v37_phase2 sdf, for the aligned pacbio reads
rule download_1000g_sdf:
    params: job_name = 'download_1000g_sdf',
    output: directory('data/genomes/1000g.sdf')
    shell:
        '''
        wget {TG_v37_SDF_URL} -O data/genomes/1000g_v37_phase2.sdf.zip;
        unzip data/genomes/1000g_v37_phase2.sdf.zip -d data/genomes
        mv data/genomes/1000g_v37_phase2.sdf data/genomes/1000g.sdf
        '''

rule download_hg38_sdf:
    params: job_name = 'download_GRCh38_sdf',
    output: directory('data/genomes/hg38.sdf')
    shell:
        '''
        wget {GRCh38_SDF_URL} -O data/genomes/GRCh38.sdf.zip;
        unzip data/genomes/GRCh38.sdf.zip -d data/genomes
        mv data/genomes/GRCh38.sdf data/genomes/hg38.sdf
        '''

rule sort_hg38_genome_track:
    params: job_name = 'sort_hg38_genome_track.{track}',
    input: track = 'genome_tracks/{track}_hg38.unsorted.bed.gz',
           genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt'
    output: track = 'genome_tracks/{track}_hg38.bed.gz',
    shell:
        '''
        {BEDTOOLS} sort -faidx {input.genome_file} -i {input.track} | \
        bgzip -c > {output.track}
        '''

rule find_coding_exons_segdups:
    params: job_name = 'find_coding_exons_segdups',
    input: exons = 'genome_tracks/coding_exons_{build}.bed.gz',
           segdup95 = 'genome_tracks/segmental_duplications_0.95_similar_{build}.bed',
           names = 'genome_tracks/{build}.chrom.sizes.txt'
    output: track = 'genome_tracks/coding_exons_intersect_segmental_duplications_{build}_0.95_similar.bed.gz'
    shell:
        '''
        bedtools intersect -a {input.exons} -b {input.segdup95} | \
        {BEDTOOLS} sort -g {input.names} | \
        {BEDTOOLS} merge | \
        bgzip -c > {output.track}
        '''

rule convert_genome_track_to_1000g:
    params: job_name = 'convert_genome_track_to_1000g.{track}',
    input: track = 'genome_tracks/{track}_hg19.bed.gz',
           names = 'genome_tracks/names.txt'
    output: track = 'genome_tracks/{track}_1000g.bed.gz'
    shell:
        '''
        gunzip -c {input.track} | \
        python3 filter_bed_chroms.py --remove_chr | \
        {BEDTOOLS} sort -g {input.names} | \
        bgzip -c > {output.track}
        '''

rule calculate_coverage:
    params: job_name = 'calculate_coverage.{individual}.{build}.{tech}.{info}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam',
            bai = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam.bai',
            genome_file = 'genome_tracks/{build}.chrom.sizes.txt',
            n_regions = 'genome_tracks/N_regions_{build}.bed',
            hg38_centromere = 'genome_tracks/sv_repeat_telomere_centromere_hg38.bed'
    output: random_pos = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam.random_pos_for_coverage',
            med_cov = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam.median_coverage',
            mean_cov = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam.mean_coverage',
            stdev_cov = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam.stdev_coverage'
    run:
        assert(wildcards.build in ['1000g','hg38'])
        # if we're using hg38, we need to remove positions in the centromeres
        centromere_subtract = '| {{BEDTOOLS}} subtract -a stdin -b {}'.format(input.hg38_centromere) if wildcards.build == 'hg38' else ''
        # sample 100000 random positions, remove N bases and centromeres (for hg38), and sort
        shell('''
        {BEDTOOLS} random -l 1 -n 100000 -g {input.genome_file} | \
        {BEDTOOLS} subtract -a stdin -b {input.n_regions} ''' + centromere_subtract + ''' |
        {BEDTOOLS} sort -faidx {input.genome_file} > {output.random_pos}''')

        add_chr = (wildcards.individual == 'NA12878' and wildcards.build == '1000g' and wildcards.tech == 'pacbio' and 'blasr' in wildcards.info)

        med_cov, mean_cov, stdev_cov = calculate_coverage(input.bam, output.random_pos, add_chr=add_chr)
        with open(output.med_cov,'w') as outf:
            print(med_cov,file=outf)
        with open(output.mean_cov,'w') as outf:
            print(mean_cov,file=outf)
        with open(output.stdev_cov,'w') as outf:
            print(stdev_cov,file=outf)

# SUBSAMPLE ILLUMINA BAM
rule subsample_illumina_60x:
    params: job_name = 'subsample_illumina_{individual}.{build}.{cov}x'
    input:  bam = 'data/{individual}.{build}/aligned_reads/illumina/illumina.aligned.all.60x.bam'
    output: bam = 'data/{individual,NA\d+}.{build}/aligned_reads/illumina/illumina.aligned.all.{cov,[0-5][0-9]}x.bam'
    run:
        subsample_frac = float(wildcards.cov) / 60.0
        shell('{SAMTOOLS} view -hb {input.bam} -s {subsample_frac} > {output.bam}')

rule tabix_index:
    params: job_name = lambda wildcards: 'tabix_index.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.{filetype}.gz'
    output: '{x}.{filetype,(bed|vcf)}.gz.tbi'
    shell:  '{TABIX} -f -p {wildcards.filetype} {input}'

# simlinks weren't working here for some reason... changed to a copy
rule make_hs37d5_alias:
    params: job_name = 'make_hs37d5_alias'
    input:  fa = 'data/genomes/hs37d5.fa',
            fai = 'data/genomes/hs37d5.fa.fai'
    output: fa = 'data/genomes/1000g.fa',
            fai = 'data/genomes/1000g.fa.fai'
    shell:  'cp {input.fa} {output.fa}; cp {input.fai} {output.fai}'

# bgzip vcf
rule bgzip_vcf_calls:
    params: job_name = 'bgzip_vcf_calls.{individual}.{build}.{calls_name}.{chrom}'
    input:  'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf'
    output: 'data/{individual}.{build}/variants/{calls_name}/{chrom,(all|\d+|X|Y)}.vcf.gz'
    shell:  '{BGZIP} -c {input} > {output}'

# bgzip vcf
rule bgzip_ground_truth:
    params: job_name = 'bgzip_ground_truth.{individual}.{build}'
    input:  'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf'
    output: 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz'
    shell:  '{BGZIP} -c {input} > {output}'

# bgzip bed
rule bgzip_bed:
    params: job_name = 'bgzip_bed.{individual}.{build}.{calls_name}.{build}'
    input: 'data/{individual}.{build}/variants/{calls_name}/region_filter.bed'
    output: 'data/{individual}.{build}/variants/{calls_name}/region_filter.bed.gz'
    shell:  '{BGZIP} -c {input} > {output}'

# gunzip fasta
rule gunzip_fasta:
    params: job_name = lambda wildcards: 'gunzip_fasta.data.genomes.{}'.format(str(wildcards.x).replace("/", "."))
    input:  'data/genomes/{x}.fa.gz'
    output: 'data/genomes/{x}.fa'
    shell:  'gunzip {input}'

# gunzip_bed
#rule gunzip_bed:
#    params: job_name = lambda wildcards: 'gunzip_bed.{id}.{build}.{}'.format(str(wildcards.x).replace("/", "."))
#    input:  'data/{id}.{build}/{x}.bed.gz'
#    output: 'data/{id}.{build,(1000g|hg38)}/{x}.bed'
#    shell:  'gunzip -c {input} > {output}'

# gunzip_bed
rule gunzip_genome_tracks_bed:
    params: job_name = 'gunzip_bed.{x}'
    input:  'genome_tracks/{x}.bed.gz'
    output: 'genome_tracks/{x}.bed'
    shell:  'gunzip -c {input} > {output}'

# index fasta reference
rule index_fasta:
    params: job_name = lambda wildcards: 'index_fa.{}'.format(str(wildcards.x).replace("/", "."))
    input:  fa  = '{x}.fa'
    output: fai = '{x}.fa.fai'
    shell:  '{SAMTOOLS} faidx {input.fa}'

# index fasta for minimap2
rule index_minimap2:
    params: job_name = lambda wildcards: 'index_minimap2.{}'.format(str(wildcards.x).replace("/", "."))
    input:  fa  = '{x}.fa'
    output: mmi = '{x}.fa.mmi'
    shell:  '{MINIMAP2} -d {output.mmi} {input.fa}'

rule index_bam:
    params: job_name = lambda wildcards: 'index_bam.{}'.format(str(wildcards.x).replace("/", "."))
    input:  bam = '{x}.bam'
    output: bai = '{x}.bam.bai'
    shell:  '{SAMTOOLS} index {input.bam} {output.bai}'

# BWA index
rule bwa_index_fasta:
    params: job_name = lambda wildcards: 'bwa_index_fasta.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.fa'
    output: '{x}.fa.bwt'
    shell: '{BWA} index {input}'

rule backup_longshot_run:
    shell:
        '''
        rm -rf data/BAK
        mkdir data/BAK

        mkdir -p data/BAK/NA12878.1000g/variants
        mkdir -p data/BAK/NA12878.1000g/vcfeval
        mv data/NA12878.1000g/variants/longshot* data/BAK/NA12878.1000g/variants 2> /dev/null
        mv data/NA12878.1000g/vcfeval/longshot* data/BAK/NA12878.1000g/vcfeval 2> /dev/null

        mkdir -p data/BAK/NA24143.hg38/variants
        mkdir -p data/BAK/NA24143.hg38/vcfeval
        mv data/NA24143.hg38/variants/longshot* data/BAK/NA24143.hg38/variants 2> /dev/null
        mv data/NA24143.hg38/vcfeval/longshot* data/BAK/NA24143.hg38/vcfeval 2> /dev/null

        mkdir -p data/BAK/NA24149.hg38/variants
        mkdir -p data/BAK/NA24149.hg38/vcfeval
        mv data/NA24149.hg38/variants/longshot* data/BAK/NA24149.hg38/variants 2> /dev/null
        mv data/NA24149.hg38/vcfeval/longshot* data/BAK/NA24149.hg38/vcfeval 2> /dev/null

        mkdir -p data/BAK/NA24385.hg38/variants
        mkdir -p data/BAK/NA24385.hg38/vcfeval
        mv data/NA24385.hg38/variants/longshot* data/BAK/NA24385.hg38/variants 2> /dev/null
        mv data/NA24385.hg38/vcfeval/longshot* data/BAK/NA24385.hg38/vcfeval 2> /dev/null

        #mkdir -p data/BAK/NA24143.1000g/variants
        #mkdir -p data/BAK/NA24143.1000g/vcfeval
        #mv data/NA24143.1000g/variants/longshot* data/BAK/NA24143.1000g/variants 2> /dev/null
        #mv data/NA24143.1000g/vcfeval/longshot* data/BAK/NA24143.1000g/vcfeval 2> /dev/null

        #mkdir -p data/BAK/NA24149.1000g/variants
        #mkdir -p data/BAK/NA24149.1000g/vcfeval
        #mv data/NA24149.1000g/variants/longshot* data/BAK/NA24149.1000g/variants 2> /dev/null
        #mv data/NA24149.1000g/vcfeval/longshot* data/BAK/NA24149.1000g/vcfeval 2> /dev/null

        #mkdir -p data/BAK/NA24385.1000g/variants
        #mkdir -p data/BAK/NA24385.1000g/vcfeval
        #mv data/NA24385.1000g/variants/longshot* data/BAK/NA24385.1000g/variants 2> /dev/null
        #mv data/NA24385.1000g/vcfeval/longshot* data/BAK/NA24385.1000g/vcfeval 2> /dev/null

        #mkdir data/BAK/plots
        #mv data/plots data/BAK/plots 2> /dev/null

        mv data/aj_trio data/BAK/aj_trio
        '''
