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
from calculate_median_coverage import calculate_median_coverage

# DATA URLs
HG19_URL     = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit'
HS37D5_URL     = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
HG38_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
HG19_SDF_URL   = 'https://s3.amazonaws.com/rtg-datasets/references/hg19.sdf.zip'
TG_v37_SDF_URL = 'https://s3.amazonaws.com/rtg-datasets/references/1000g_v37_phase2.sdf.zip'
GRCh38_SDF_URL = 'https://s3.amazonaws.com/rtg-datasets/references/GRCh38.sdf.zip'

# PATHS TO TOOLS
FASTQUTILS = '/home/pedge/installed/ngsutils/bin/fastqutils'
TWOBITTOFASTA = 'twoBitToFa' # can be downloaded from 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'
SAMTOOLS       = 'samtools' # v1.3
FASTQ_DUMP     = 'fastq-dump' # v2.5.2
REAPER         = '../target/release/reaper' # v0.1
RTGTOOLS       = 'rtg' #'/home/pedge/installed/rtg-tools-3.8.4/rtg' # v3.8.4, https://www.realtimegenomics.com/products/rtg-tools
BGZIP = 'bgzip'
TABIX = 'tabix'
FREEBAYES      = '/home/pedge/git/freebayes/bin/freebayes'
SIMLORD = 'simlord'
DWGSIM = '/home/pedge/git/DWGSIM/dwgsim'
BWA = '/home/pedge/installed/bwa'
BLASR = 'blasr'
BAX2BAM = 'bax2bam'
SAWRITER = 'sawriter'
NGMLR = 'ngmlr'
MINIMAP2 = 'minimap2'
BCFTOOLS = '/opt/biotools/bcftools/bin/bcftools'
PYFAIDX = '/home/pedge/installed/opt/python/bin/faidx'
BEDTOOLS = 'bedtools' # v 2.27
EXTRACTHAIRS = 'HapCUT2/build/extractHAIRS'
HAPCUT2 = 'HapCUT2/build/HAPCUT2'
PYTHON = 'python3'
MAP_COUNTER = 'map_counter/target/release/map_counter'

chroms = ['{}'.format(i) for i in range(1,23)]
ref_file = {'1000g':'data/genomes/hs37d5.fa', 'hg38':'data/genomes/hg38.fa'}

include: "simulation.1000g.snakefile"
include: "NA12878.1000g.snakefile"
include: "NA24385.hg38.snakefile"  # AJ Son,    hg38
include: "NA24143.hg38.snakefile"  # AJ Mother, hg38
include: "NA24149.hg38.snakefile"  # AJ Father, hg38
include: "NA24385.1000g.snakefile" # AJ Son,    1000g
include: "NA24143.1000g.snakefile" # AJ Mother, 1000g
include: "NA24149.1000g.snakefile" # AJ Father, 1000g
include: "aj_trio.snakefile" #
include: "paper_tables_and_figures.snakefile"
include: "haplotyping.snakefile"
include: "make_variant_counts_table.snakefile"
include: "duplicated_gene_visualization.snakefile"

# DEFAULT
rule all:
    input:
        'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.median_coverage',
        #'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam.median_coverage',
        expand('data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.{cov}x.bam.median_coverage',cov=[20,30,40,50,69]),
        'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
        'data/NA24149.hg38/aligned_reads/pacbio/pacbio.blasr.all.32x.bam.median_coverage',
        #'data/plots/NA24149.hg38_prec_recall_all.png',
        #'data/plots/NA24143.hg38_prec_recall_all.png',
        #'data/plots/NA24385.hg38_prec_recall_all.png',
        #'data/plots/NA12878.1000g_prec_recall_all.png',
        #'data/plots/plot_mappability_bars.simulation.1000g.png'
        #expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac{mapfrac}/pacbio.blasr.60x__all.map_count.txt',cov=[10,20,30,40,50],mapfrac=[0.0,0.5,0.75,0.9,1.0]),
        #expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac{mapfrac}/illumina.60x__all.map_count.txt',cov=[10,20,30,40,50],mapfrac=[0.0,0.5,0.75,0.9,1.0])
        #'data/simulation.1000g/aligned_reads/pacbio/segdup95_stats/all.stats.bed'
        #'data/plots/plot_fp_near_indel.NA24385.hg38.png'
        #expand('data/output/variant_analysis_fp_fn__NA24385.hg38__GQ50__reaper.pacbio.blasr.{cov}x.-z__1.tex',cov=[20,30,40,50,69]),
        #'data/output/variant_analysis_fp_fn__NA12878.1000g__GQ44__reaper.pacbio.blasr.44x.-z__1.tex'
        #'data/output/variant_counts_table.NA12878.1000g.blasr.44.GQ44.tex'
        #'data/NA12878.1000g/duplicated_gene_vis/pacbio.44x.dup_genes.bam.bai',   # bam file for duplicated gene visualization
        #'data/NA12878.1000g/duplicated_gene_vis/illumina.30x.dup_genes.bam.bai'  # bam file for duplicated gene visualization
        #'data/output/four_GIAB_genomes_table.aj_trio_hg38_blasr.all.tex',
        #'data/output/four_GIAB_genomes_table_extended.aj_trio_hg38_blasr.all.tex',
        #'data/plots/haplotyping_results_barplot.png'
        #'data/plots/fig3_precision_recall_bars_NA24385_NA12878.blasr.hg38.png'
        #'data/plots/depth_vs_breadth_mappability.simulation.30x.png',
        #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/NA24143.30x.pileup_over_mendelian_positions.txt',
        #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/NA24149.32x.pileup_over_mendelian_positions.txt',
        #'data/plots/simulation_pr_barplot_genome_vs_segdup.all.GQ50.png',
        #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/all.vcf.gz',
        #'data/output/variant_counts_table.NA12878.1000g.blasr.44.GQ44.tex'

        #'data/plots/NA12878.1000g_prec_recall_all.png',
        #'data/plots/NA24143.hg38_prec_recall_all.png',
        #'data/plots/NA24149.hg38_prec_recall_all.png',
        #'data/plots/NA24385.hg38_prec_recall_all.png',
        #'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/mendelian/all.vcf.gz',
        #'data/plots/simulation_pr_barplot_genome_vs_segdup.all.GQ50.png',
        #'data/plots/simulation_pr_barplot_genome_vs_segdup_extended.1.GQ50.png',
        #'data/plots/depth_vs_breadth_mappability.NA12878.30x.png',
        #'data/plots/NA24385.1000g_prec_recall_all.png',
        #'data/plots/NA24143.1000g_prec_recall_all.png',
        #'data/plots/NA24149.1000g_prec_recall_all.png',
        #'data/plots/supp_fig3_with_haplotypefree_precision_recall_bars_NA24385_NA12878.bwamem.1000g.png',
        #'data/plots/supp_fig3_with_haplotypefree_precision_recall_bars_NA24385_NA12878.blasr.hg38.png',

        #'data/output/four_GIAB_genomes_table.aj_trio_1000g_bwamem.all.GQ50.tex', # current version with bwamem
        #'data/output/four_GIAB_genomes_table_extended.aj_trio_1000g_bwamem.all.GQ50.tex', # current version with bwamem

        #'data/plots/effect_of_haplotyping.giab_individuals.prec_recall_all.png',
        #'data/output/four_GIAB_genomes_table.aj_trio_hg38_blasr.all.GQ50.tex', # once hg38 blasr bams available for aj mother + father
        #'data/output/four_GIAB_genomes_table_extended.aj_trio_hg38_blasr.all.GQ50.tex', # once hg38 blasr bams available for aj mother + father
        #'data/plots/simulation.1_prec_recall_1.png',
        #'data/plots/simulation_prec_recall_in_segdups_1.png',
        #'data/plots/chr1_simulated_60x_pacbio_mismapped_read_distribution.segdup.png',
        #'data/plots/chr1_simulated_60x_pacbio_mismapped_read_distribution.png',
        #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/1.vcf.gz',
        #'data/output/four_GIAB_genomes_table.1.GQ50.tex',
        #'data/plots/simulation_pr_barplot_genome_vs_segdup.1.GQ50.png',
        #'data/plots/effect_of_haplotyping.giab_individuals.prec_recall_1.png',
        #'data/plots/PR_curve_3_mappers_AJ_father_chr20.png',
        #'data/plots/compare_mappers_reaper_in_segdups_simulation_1.png'
        #'data/plots/compare_mappers_reaper_in_segdups_simulation_1.png',
        #'data/plots/simulation_prec_recall_ngmlr_1.png',
        #'data/output/variant_analysis_fp_fn__NA12878__reaper.pacbio.blasr.44x.-z__1.tex',
        #'data/output/variant_analysis_fp_fn__NA24385__reaper.pacbio.ngmlr.69x.-z__1.tex',
        #'data/output/variant_analysis_fp_fn__NA24149__reaper.pacbio.ngmlr.32x.-z__1.tex',
        #'data/output/variant_analysis_fp_fn__NA24143__reaper.pacbio.ngmlr.30x.-z__1.tex',
        #'data/plots/NA12878_prec_recall_all.png',
        #'data/NA12878.1000g/vcfeval/illumina_30x.filtered/all.done',
        #'data/NA24385.1000g/vcfeval/illumina_30x.filtered/all.done',
        #'data/NA24143.1000g/vcfeval/illumina_30x.filtered/all.done',
        #'data/NA24149.1000g/vcfeval/illumina_30x.filtered/all.done',
        #'data/NA24385.hg38/vcfeval/illumina_30x.filtered/all.done',
        #'data/NA24143.hg38/vcfeval/illumina_30x.filtered/all.done',
        #'data/NA24149.hg38/vcfeval/illumina_30x.filtered/all.done',


rule run_pileups:
    params: job_name = 'run_pileups_{individual}.{cov}x',
    input:  bam = 'data/{individual}.hg38/aligned_reads/pacbio/pacbio.blasr.all.{cov}x.bam',
            pos = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/all.positions_only.txt'
    output: txt = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/{individual}.{cov}x.pileup_over_mendelian_positions.txt',
    shell: 'python3 run_pileups.py {input.bam} {input.pos} > {output.txt}'

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

rule rtg_filter_reaper_MEC_AQ:
    params: job_name = 'rtg_filter_reaper_MEC_AQ.{individual}.{build}.{aligner}.{cov}x.chr{chrom}',
    input:  vcfgz = 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{cov}x.-z/{chrom}.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{cov}x.-z/{chrom}.vcf.gz.tbi',
    output: vcf = 'data/{individual}.{build}/variants/reaper.filtered.pacbio.{aligner}.{cov}x.-z/{chrom}.vcf',
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter --keep-expr "INFO.AQ > 7" --fail=AQ -i {input.vcfgz} -o {output.vcf}.gz
        gunzip {output.vcf}.gz
        '''

rule vcfeval_rtgtools:
    params: job_name = 'vcfeval_rtgtools.{individual}.{build}.{calls_name}.{chrom}',
            region_arg = lambda wildcards: '--region={}'.format(chr_prefix(wildcards.chrom,wildcards.build)) if wildcards.chrom != 'all' else ''
    input:  calls_vcf = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf.gz',
            calls_ix = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf.gz.tbi',
            ground_truth = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz',
            ground_truth_ix = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz.tbi',
            region_filter ='data/{individual}.{build}/variants/ground_truth/region_filter.bed',
            sdf = 'data/genomes/{build}.sdf'
    output: done = 'data/{individual}.{build}/vcfeval/{calls_name}/{chrom}.done'
    shell:
        '''
        rm -rf data/{wildcards.individual}.{wildcards.build}/vcfeval/{wildcards.calls_name}/{wildcards.chrom}
        {RTGTOOLS} RTG_MEM=12g vcfeval \
        {params.region_arg} \
        -c {input.calls_vcf} \
        -b {input.ground_truth} \
        -e {input.region_filter} \
        -t {input.sdf} \
        -o data/{wildcards.individual}.{wildcards.build}/vcfeval/{wildcards.calls_name}/{wildcards.chrom};
        cp data/{wildcards.individual}.{wildcards.build}/vcfeval/{wildcards.calls_name}/{wildcards.chrom}/done {output.done};
        '''

rule rtg_decompose_variants_ground_truth:
    params: job_name = 'rtg_decompose_variants_ground_truth.{individual}.{build}',
    input:  vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz.tbi'
    output: vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfdecompose --break-mnps --break-indels -i {input.vcfgz} -o {output.vcfgz}'

rule analyze_variants:
    params: job_name = 'analyze_variants.{individual}.{build}.{chrom}.{calls_name}.GQ{GQ}',
    input:  vcfeval = 'data/{individual}.{build}/vcfeval/{calls_name}/{chrom}.done',
            pacbio_calls_vcfgz = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.GQ{GQ}.vcf.gz',
            pacbio_calls_ix = 'data/{individual}.{build}/variants/{calls_name}/{chrom}.GQ{GQ}.vcf.gz.tbi',
            ground_truth_vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz',
            ground_truth_ix = 'data/{individual}.{build}/variants/ground_truth/ground_truth.vcf.gz.tbi',
            ground_truth_bed = 'data/{individual}.{build}/variants/ground_truth/region_filter.bed.gz',
            ground_truth_bed_ix = 'data/{individual}.{build}/variants/ground_truth/region_filter.bed.gz.tbi',
            random_positions_vcfgz = 'data/{individual}.{build}/random_positions/{chrom}.confident_only.vcf.gz',
            random_positions_ix = 'data/{individual}.{build}/random_positions/{chrom}.confident_only.vcf.gz.tbi',
            str_bed = 'genome_tracks/STRs_{build}.bed.gz',
            str_bed_ix = 'genome_tracks/STRs_{build}.bed.gz.tbi',
            line_bed = 'genome_tracks/LINEs_{build}.bed.gz',
            line_bed_ix = 'genome_tracks/LINEs_{build}.bed.gz.tbi',
            sine_bed = 'genome_tracks/SINEs_{build}.bed.gz',
            sine_bed_ix = 'genome_tracks/SINEs_{build}.bed.gz.tbi',
            ref_1000g_fa = 'data/genomes/hs37d5.fa',
            ref_1000g_fai = 'data/genomes/hs37d5.fa.fai',
            ref_hg38_fa = 'data/genomes/hg38.fa',
            ref_hg38_fai = 'data/genomes/hg38.fa.fai'
    output: tex = 'data/output/variant_analysis_fp_fn__{individual}.{build}__GQ{GQ}__{calls_name}__{chrom}.tex'
    run:
        analyze_variants(chrom_name = chr_prefix(wildcards.chrom, wildcards.build),
                         pacbio_calls_vcfgz = input.pacbio_calls_vcfgz,
                         fp_calls_vcfgz = os.path.join(input.vcfeval[:-5],'fp.vcf.gz'),
                         fn_calls_vcfgz =  os.path.join(input.vcfeval[:-5],'fn.vcf.gz'),
                         ground_truth_vcfgz = input.ground_truth_vcfgz,
                         ground_truth_bed_file = input.ground_truth_bed,
                         random_positions_vcfgz = input.random_positions_vcfgz,
                         str_tabix_bed_file = input.str_bed,
                         line_tabix_bed_file = input.line_bed,
                         sine_tabix_bed_file = input.sine_bed,
                         ref_fa = ref_file[wildcards.build],
                         gq_cutoff = float(wildcards.GQ),
                         output_file = output.tex)

def get_chr_len(chrom, sizes_file):
    with open(sizes_file,'r') as inf:
        for line in inf:
            el = line.strip().split('\t')
            if chrom == el[0]:
                return int(el[1])

        print("CHROMOSOME SIZE NOT FOUND")
        exit(1)

rule filter_random_positions:
    params: job_name = 'filter_random_positions.{build}.{chrom}',
    input:  vcfgz = 'data/{individual}.{build}/random_positions/{chrom}.whole_chrom.vcf.gz',
            tbi   = 'data/{individual}.{build}/random_positions/{chrom}.whole_chrom.vcf.gz.tbi',
            bed   = 'data/{individual}.{build}/variants/ground_truth/region_filter.bed'
    output: vcfgz = 'data/{individual}.{build,(1000g|hg38)}/random_positions/{chrom}.confident_only.vcf.gz'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.bed} -i {input.vcfgz} -o {output.vcfgz}'

rule generate_random_positions:
    params: job_name = 'generate_random_positions.{chrom}.{build}'
    input: sizes_file = 'genome_tracks/{build}.chrom.sizes.txt'
    output: vcf = 'data/{individual}.{build,(1000g|hg38)}/random_positions/{chrom}.whole_chrom.vcf',
            vcfgz = 'data/{individual}.{build,(1000g|hg38)}/random_positions/{chrom}.whole_chrom.vcf.gz',
    run:
        chrom_name = chr_prefix(wildcards.chrom, wildcards.build)
        chrlen = get_chr_len(chrom_name, input.sizes_file)
        generate_random_calls(chrom_name, chrlen, 100000, output.vcf)
        shell('bgzip -c {output.vcf} > {output.vcfgz}')

rule rtg_filter_SNVs_ground_truth:
    params: job_name = 'rtg_filter_SNVs_ground_truth.{individual}.{build}',
    input:  vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.vcf.gz.tbi'
    output: vcfgz = 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --snps-only -i {input.vcfgz} -o {output.vcfgz}'

from filter_SNVs import filter_SNVs
rule filter_illumina_SNVs:
    params: job_name = 'filter_SNVs_illumina.{individual}.{build}.chr{chrom}',
    input:  vcf = 'data/{individual}.{build}/variants/illumina_{cov}x/{chrom}.vcf'
    output: vcf = 'data/{individual}.{build}/variants/illumina_{cov}x.filtered/{chrom,(\d+)}.vcf'
    run:
        cov_filter = int(float(wildcards.cov)*2)
        filter_SNVs(input.vcf, output.vcf, cov_filter, density_count=10, density_len=500, density_qual=50)

rule combine_chrom:
    params: job_name = 'combine_chroms.{individual}.{build}.{calls_name}',
    input: expand('data/{{individual}}.{{build}}/variants/{{calls_name}}/{chrom}.vcf',chrom=chroms)
    output: 'data/{individual}.{build}/variants/{calls_name,(illumina|reaper.pacbio).+}/all.vcf'
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

rule run_reaper:
    params: job_name = 'reaper.pacbio.{aligner}.{individual}.{build}.cov{cov}.{options}.chr{chrom}',
    input:  bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.{chrom}.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.{chrom}.{cov}x.bam.bai',
            hg19    = 'data/genomes/hg19.fa',
            hg19_ix = 'data/genomes/hg19.fa.fai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+)}.vcf',
            #debug = 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom}.debug',
            #no_hap_vcf = 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom}.debug/2.0.realigned_genotypes.vcf',
            runtime = 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{cov,\d+}x.{options}/{chrom,(\d+)}.vcf.runtime'
    run:
        options_str = wildcards.options.replace('_',' ')
        if wildcards.individual == 'NA12878':
            t1 = time.time()
            shell('{REAPER} -r chr{wildcards.chrom} -F -d {output.debug} {options_str} -s {wildcards.individual} --bam {input.bam} --ref {input.hg19} --out {output.vcf}.tmp')
            t2 = time.time()
            # remove 'chr' from reference name in vcf
            remove_chr_from_vcf(output.vcf+'.tmp',output.vcf)
            # remove 'chr' from no-haplotype version of vcf
            remove_chr_from_vcf(output.no_hap_vcf, output.no_hap_vcf+'.tmp')
            shell('mv {output.no_hap_vcf}.tmp {output.no_hap_vcf}')
        else:
            w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
            w_ref = ref_file[wildcards.build]
            t1 = time.time()
            shell('{REAPER} -r {w_chrom} -F -d {output.debug} {options_str} -s {wildcards.individual} --bam {input.bam} --ref {w_ref} --out {output.vcf}')
            t2 = time.time()

        with open(output.runtime,'w') as outf:
            print(seconds_to_formatted_time(t2-t1),file=outf)

# Call 30x Illumina variants
rule call_variants_Illumina:
    params: job_name = 'call_illumina.{individual}.{build}.{cov}x.chr{chrom}',
    input: bam = 'data/{individual}.{build}/aligned_reads/illumina/illumina.{cov}x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/illumina/illumina.{cov}x.bam.bai',
            ref_1000g_fa = 'data/genomes/hs37d5.fa',
            ref_1000g_fai = 'data/genomes/hs37d5.fa.fai',
            ref_hg38_fa = 'data/genomes/hg38.fa',
            ref_hg38_fai = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/variants/illumina_{cov}x/{chrom,(\d+)}.vcf',
            runtime = 'data/{individual}.{build}/variants/illumina_{cov}x/{chrom,(\d+)}.vcf.runtime'
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
    input:  expand('data/{{individual}}.{{build}}/aligned_reads/{{tech}}/genomecov_histograms_mapq{{mapq}}/{{info}}__{chrom}.txt', chrom=chroms)
    output: 'data/{individual}.{build}/aligned_reads/{tech}/genomecov_histograms_mapq{mapq}/{info}__all.txt'
    shell: 'cat {input} > {output}'

rule generate_coverage_histogram_illumina:
    params: job_name = 'generate_coverage_histogram_illumina.{individual}.{build}.MAPQ{mapq}.{chrom}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/illumina/illumina.30x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/illumina/illumina.30x.bam.bai'
    output: 'data/{individual}.{build}/aligned_reads/illumina/genomecov_histograms_mapq{mapq}/illumina.30x__{chrom,(\d+)}.txt'
    run:
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        shell('''
        {SAMTOOLS} view -F 3844 -q {wildcards.mapq} {input.bam} -hb {w_chrom} | \
        {BEDTOOLS} genomecov -ibam - | \
        grep -P '^{w_chrom}\\t' > {output}
        ''')

rule generate_coverage_histogram_pacbio:
    params: job_name = 'generate_coverage_histogram_pacbio.{individual}.{build}.MAPQ{mapq}.{chrom}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.blasr.{chrom}.30x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.blasr.{chrom}.30x.bam.bai'
    output: 'data/{individual}.{build}/aligned_reads/pacbio/genomecov_histograms_mapq{mapq}/pacbio.blasr.all.30x__{chrom,(\d+)}.txt'
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

rule download_HS37D5:
    params: job_name = 'download_hs37d5',
    output: 'data/genomes/hs37d5.fa'
    shell: 'wget {HS37D5_URL} -O {output}.gz; gunzip {output}.gz'

# download hg19 reference, for the aligned pacbio reads
rule download_hg19_sdf:
    params: job_name = 'download_hg19_sdf',
    output: 'data/genomes/hg19.sdf'
    shell:
        '''
        wget {HG19_SDF_URL} -O {output}.zip;
        unzip {output}.zip -d data/genomes
        '''

# download 1000g_v37_phase2 sdf, for the aligned pacbio reads
rule download_1000g_sdf:
    params: job_name = 'download_1000g_sdf',
    output: 'data/genomes/1000g.sdf'
    shell:
        '''
        wget {TG_v37_SDF_URL} -O data/genomes/1000g_v37_phase2.sdf.zip;
        unzip data/genomes/1000g_v37_phase2.sdf.zip -d data/genomes
        mv data/genomes/1000g_v37_phase2.sdf data/genomes/1000g.sdf
        '''

rule download_hg38_sdf:
    params: job_name = 'download_GRCh38_sdf',
    output: 'data/genomes/hg38.sdf'
    shell:
        '''
        wget {GRCh38_SDF_URL} -O data/genomes/GRCh38.sdf.zip;
        unzip data/genomes/GRCh38.sdf.zip -d data/genomes
        mv data/genomes/GRCh38.sdf data/genomes/hg38.sdf
        '''

rule sort_hg38_genome_track:
    params: job_name = 'sort_hg38_genome_track.{track}',
    input: track = 'genome_tracks/{track}_hg38.unsorted.bed.gz',
           genome_file = 'genome_tracks/hg38.chrom.sizes.txt'
    output: track = 'genome_tracks/{track}_hg38.bed.gz',
    shell:
        '''
        {BEDTOOLS} sort -faidx {input.genome_file} -i {input.track} | \
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

# count_bam_median_coverage
rule calculate_median_coverage:
    params: job_name = 'calculate_median_coverage.{individual}.{build}.{tech}.{info}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam',
            bai = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam',
            genome_file = 'genome_tracks/{build}.chrom.sizes.txt'
    output: random_pos = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam.for_median_coverage.random_pos',
            cov = 'data/{individual}.{build}/aligned_reads/{tech}/{info}.bam.median_coverage'
    run:
        shell('{BEDTOOLS} random -l 1 -n 10000 -g {input.genome_file} > {output.random_pos}')
        add_chr = (wildcards.individual == 'NA12878' and wildcards.build == '1000g' and wildcards.tech == 'pacbio')

        med_cov = calculate_median_coverage(input.bam, output.random_pos, add_chr=add_chr)
        with open(output.cov,'w') as outf:
            print(med_cov,file=outf)

# SUBSAMPLE ILLUMINA BAM
rule subsample_illumina_60x:
    params: job_name = 'subsample_illumina_{individual}.{build}.{cov}x'
    input:  'data/{individual}.{build}/aligned_reads/illumina/illumina.60x.bam'
    output: 'data/{individual,NA\d+}.{build}/aligned_reads/illumina/illumina.{cov,(1|2|3|4|5)0}x.bam'
    run:
        subsample_frac = float(wildcards.cov) / 60.0
        shell('{SAMTOOLS} view -hb {input} -s {subsample_frac} > {output}')

rule tabix_index:
    params: job_name = lambda wildcards: 'tabix_index.{}'.format(str(wildcards.x).replace("/", "."))
    input:  '{x}.{filetype}.gz'
    output: '{x}.{filetype,(bed|vcf)}.gz.tbi'
    shell:  '{TABIX} -f -p {wildcards.filetype} {input}'

# bgzip vcf
rule bgzip_vcf_calls:
    params: job_name = 'bgzip_vcf_calls.{individual}.{build}.{calls_name}.{chrom}'
    input:  'data/{individual}.{build}/variants/{calls_name}/{chrom}.vcf'
    output: 'data/{individual}.{build}/variants/{calls_name}/{chrom,(all|\d+)}.vcf.gz'
    shell:  '{BGZIP} -c {input} > {output}'

# bgzip vcf
rule bgzip_debug_vcf_calls:
    params: job_name = 'bgzip_no_hap_vcf_calls.{x}.{individual}.{build}.{calls_name}.{chrom}'
    input:  'data/{individual}.{build}/variants/{calls_name}/{chrom}.debug/{x}.vcf'
    output: 'data/{individual}.{build}/variants/{calls_name}/{chrom,(all|\d+)}.debug/{x}.vcf.gz'
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

rule backup_reaper_run:
    shell:
        '''
        rm -rf data/BAK
        mkdir data/BAK

        mkdir -p data/BAK/NA12878.1000g/variants
        mkdir -p data/BAK/NA12878.1000g/vcfeval
        mv data/NA12878.1000g/variants/reaper* data/BAK/NA12878.1000g/variants 2> /dev/null
        mv data/NA12878.1000g/vcfeval/reaper* data/BAK/NA12878.1000g/vcfeval 2> /dev/null

        mkdir -p data/BAK/NA24143.hg38/variants
        mkdir -p data/BAK/NA24143.hg38/vcfeval
        mv data/NA24143.hg38/variants/reaper* data/BAK/NA24143.hg38/variants 2> /dev/null
        mv data/NA24143.hg38/vcfeval/reaper* data/BAK/NA24143.hg38/vcfeval 2> /dev/null

        mkdir -p data/BAK/NA24149.hg38/variants
        mkdir -p data/BAK/NA24149.hg38/vcfeval
        mv data/NA24149.hg38/variants/reaper* data/BAK/NA24149.hg38/variants 2> /dev/null
        mv data/NA24149.hg38/vcfeval/reaper* data/BAK/NA24149.hg38/vcfeval 2> /dev/null

        mkdir -p data/BAK/NA24385.hg38/variants
        mkdir -p data/BAK/NA24385.hg38/vcfeval
        mv data/NA24385.hg38/variants/reaper* data/BAK/NA24385.hg38/variants 2> /dev/null
        mv data/NA24385.hg38/vcfeval/reaper* data/BAK/NA24385.hg38/vcfeval 2> /dev/null

        #mkdir -p data/BAK/NA24143.1000g/variants
        #mkdir -p data/BAK/NA24143.1000g/vcfeval
        #mv data/NA24143.1000g/variants/reaper* data/BAK/NA24143.1000g/variants 2> /dev/null
        #mv data/NA24143.1000g/vcfeval/reaper* data/BAK/NA24143.1000g/vcfeval 2> /dev/null

        #mkdir -p data/BAK/NA24149.1000g/variants
        #mkdir -p data/BAK/NA24149.1000g/vcfeval
        #mv data/NA24149.1000g/variants/reaper* data/BAK/NA24149.1000g/variants 2> /dev/null
        #mv data/NA24149.1000g/vcfeval/reaper* data/BAK/NA24149.1000g/vcfeval 2> /dev/null

        #mkdir -p data/BAK/NA24385.1000g/variants
        #mkdir -p data/BAK/NA24385.1000g/vcfeval
        #mv data/NA24385.1000g/variants/reaper* data/BAK/NA24385.1000g/variants 2> /dev/null
        #mv data/NA24385.1000g/vcfeval/reaper* data/BAK/NA24385.1000g/vcfeval 2> /dev/null

        #mkdir data/BAK/plots
        #mv data/plots data/BAK/plots 2> /dev/null

        mv data/aj_trio data/BAK/aj_trio
        '''
