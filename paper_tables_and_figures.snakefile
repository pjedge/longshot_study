from paper_tables_and_figures import genomes_table_files
from matplotlib_venn import venn3

#rule plot_true_indel_FP_fixed_GQ:
#    params: job_name = 'plot_true_indel_FP_fixed_GQ',
#    input:
#        expand('data/output/variant_analysis_fp_fn__NA24385.hg38__GQ50__longshot.pacbio.blasr.{cov}x.___1.tex',cov=[20,30,40,50,69])

rule plot_pr_curve_4panel:
    params: job_name = 'plot_pr_curve_4panel.{aligner}.{chrom}'
    input:
        NA12878_il30 = 'data/NA12878.1000g/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
        NA12878_pb30 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.30x._/{chrom}',
        NA12878_pb44 = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.44x._/{chrom}',
        NA12878_il30_cov = 'data/NA12878.1000g/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
        NA12878_pb30_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage',
        NA12878_pb44_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.{aligner}.all.44x.bam.median_coverage',
        NA24385_il30 = 'data/NA24385.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
        NA24385_pb30 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.30x._/{chrom}',
        NA24385_pb40 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.40x._/{chrom}',
        NA24385_pb50 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.50x._/{chrom}',
        NA24385_pb69 = 'data/NA24385.hg38/vcfeval/longshot.pacbio.{aligner}.69x._/{chrom}',
        NA24385_il30_cov = 'data/NA24385.hg38/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
        NA24385_pb30_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage',
        NA24385_pb40_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.40x.bam.median_coverage',
        NA24385_pb50_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.50x.bam.median_coverage',
        NA24385_pb69_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.69x.bam.median_coverage',
        NA24149_il34 = 'data/NA24149.hg38/vcfeval/freebayes.illumina.aligned.34x.filtered/{chrom}',
        NA24149_pb32 = 'data/NA24149.hg38/vcfeval/longshot.pacbio.{aligner}.32x._/{chrom}',
        NA24149_il34_cov = 'data/NA24149.hg38/aligned_reads/illumina/illumina.aligned.all.34x.bam.median_coverage',
        NA24149_pb32_cov = 'data/NA24149.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.32x.bam.median_coverage',
        NA24143_il30 = 'data/NA24143.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/{chrom}',
        NA24143_pb30 = 'data/NA24143.hg38/vcfeval/longshot.pacbio.{aligner}.30x._/{chrom}',
        NA24143_il30_cov = 'data/NA24143.hg38/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
        NA24143_pb30_cov = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.{aligner}.all.30x.bam.median_coverage'
    output:
        png = 'data/plots/prec_recall_4panel_{aligner}.{chrom}.png'
    run:
        ptf.plot_vcfeval_4panel(dirlist1=[input.NA12878_il30, input.NA12878_pb30, input.NA12878_pb44],
                          labels1=['Freebayes, Illumina {}'.format(parse_int_file(input.NA12878_il30_cov)) + r'$\times$',
                          'Longshot, PacBio {}'.format(parse_int_file(input.NA12878_pb30_cov)) + r'$\times$',
                          'Longshot, PacBio {}'.format(parse_int_file(input.NA12878_pb44_cov)) + r'$\times$'],
                          colors1=['r','#8080ff','#3333ff'],
                          legendloc1='lower left',
                          title1='NA12878',
                          dirlist2=[input.NA24385_il30, input.NA24385_pb30, input.NA24385_pb40, input.NA24385_pb50, input.NA24385_pb69],
                          labels2=['Freebayes, Illumina {}'.format(parse_int_file(input.NA24385_il30_cov)) + r'$\times$',
                          'Longshot, PacBio {}'.format(parse_int_file(input.NA24385_pb30_cov)) + r'$\times$',
                          'Longshot, PacBio {}'.format(parse_int_file(input.NA24385_pb40_cov)) + r'$\times$',
                          'Longshot, PacBio {}'.format(parse_int_file(input.NA24385_pb50_cov)) + r'$\times$',
                          'Longshot, PacBio {}'.format(parse_int_file(input.NA24385_pb69_cov)) + r'$\times$'],
                          colors2=['#ff0707','#8080ff','#6666ff','#3333ff','b'],
                          legendloc2='lower left',
                          title2='NA24385',
                          dirlist3=[input.NA24149_il34, input.NA24149_pb32],
                          labels3=['Freebayes, Illumina {}'.format(parse_int_file(input.NA24149_il34_cov)) + r'$\times$',
                          'Longshot, PacBio {}'.format(parse_int_file(input.NA24149_pb32_cov)) + r'$\times$'],
                          colors3=['r','#3333ff','#ccccff','#9999ff','#8080ff','#6666ff'],
                          legendloc3='lower left',
                           title3='NA24149',
                           dirlist4=[input.NA24143_il30, input.NA24143_pb30],
                          labels4=['Freebayes, Illumina {}'.format(parse_int_file(input.NA24143_il30_cov)) + r'$\times$',
                                   'Longshot, PacBio {}'.format(parse_int_file(input.NA24143_pb30_cov)) + r'$\times$'] ,
                           colors4=['r','#3333ff','#ccccff','#9999ff','#8080ff','#6666ff'],
                          legendloc4='lower left',
                           title4='NA24143',
                           xlim=(0.6,1.0),
                           ylim=(0.98599,1.0),
                           output_file=output.png
                          )

rule make_four_genomes_table_extended:
    params: job_name = 'make_four_genomes_table_extended.aj_trio_{aj_trio_build}_{aj_trio_aligner}.{chrom}'
    input:
        NA12878_30x_vcfeval = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.30x._/{chrom}',
        NA12878_30x_vcfstats_genome = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.30x._/{chrom}.GQ30.vcf.stats',
        NA12878_30x_vcfstats_outside_giab = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.30x._/{chrom}.outside_GIAB.GQ30.vcf.stats',
        NA12878_30x_runtime = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.30x._/{chrom}.vcf.runtime',
        NA12878_30_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
        NA12878_44x_vcfeval = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.44x._/{chrom}',
        NA12878_44x_vcfstats_genome = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/{chrom}.GQ44.vcf.stats',
        NA12878_44x_vcfstats_outside_giab = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/{chrom}.outside_GIAB.GQ44.vcf.stats',
        NA12878_44x_runtime = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/{chrom}.vcf.runtime',
        NA12878_44_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.median_coverage',
        NA24385_30x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}',
        NA24385_30x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}.GQ30.vcf.stats',
        NA24385_30x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}.outside_GIAB.GQ30.vcf.stats',
        NA24385_30x_runtime = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}.vcf.runtime',
        NA24385_30_cov = 'data/NA24385.{aj_trio_build}/aligned_reads/pacbio/pacbio.{aj_trio_aligner}.all.30x.bam.median_coverage',
        NA24385_40x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/longshot.pacbio.{aj_trio_aligner}.40x._/{chrom}',
        NA24385_40x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.40x._/{chrom}.GQ40.vcf.stats',
        NA24385_40x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.40x._/{chrom}.outside_GIAB.GQ40.vcf.stats',
        NA24385_40x_runtime = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.40x._/{chrom}.vcf.runtime',
        NA24385_40_cov = 'data/NA24385.{aj_trio_build}/aligned_reads/pacbio/pacbio.{aj_trio_aligner}.all.40x.bam.median_coverage',
        NA24385_50x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/longshot.pacbio.{aj_trio_aligner}.50x._/{chrom}',
        NA24385_50x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.50x._/{chrom}.GQ50.vcf.stats',
        NA24385_50x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.50x._/{chrom}.outside_GIAB.GQ50.vcf.stats',
        NA24385_50x_runtime = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.50x._/{chrom}.vcf.runtime',
        NA24385_50_cov = 'data/NA24385.{aj_trio_build}/aligned_reads/pacbio/pacbio.{aj_trio_aligner}.all.50x.bam.median_coverage',
        NA24385_69x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/longshot.pacbio.{aj_trio_aligner}.69x._/{chrom}',
        NA24385_69x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.69x._/{chrom}.GQ69.vcf.stats',
        NA24385_69x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.69x._/{chrom}.outside_GIAB.GQ69.vcf.stats',
        NA24385_69x_runtime = 'data/NA24385.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.69x._/{chrom}.vcf.runtime',
        NA24385_69_cov = 'data/NA24385.{aj_trio_build}/aligned_reads/pacbio/pacbio.{aj_trio_aligner}.all.69x.bam.median_coverage',
        NA24149_32x_vcfeval = 'data/NA24149.{aj_trio_build}/vcfeval/longshot.pacbio.{aj_trio_aligner}.32x._/{chrom}',
        NA24149_32x_vcfstats_genome = 'data/NA24149.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.32x._/{chrom}.GQ32.vcf.stats',
        NA24149_32x_vcfstats_outside_giab = 'data/NA24149.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.32x._/{chrom}.outside_GIAB.GQ32.vcf.stats',
        NA24149_32x_runtime = 'data/NA24149.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.32x._/{chrom}.vcf.runtime',
        NA24149_32_cov = 'data/NA24149.{aj_trio_build}/aligned_reads/pacbio/pacbio.{aj_trio_aligner}.all.32x.bam.median_coverage',
        NA24143_30x_vcfeval = 'data/NA24143.{aj_trio_build}/vcfeval/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}',
        NA24143_30x_vcfstats_genome = 'data/NA24143.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}.GQ30.vcf.stats',
        NA24143_30x_vcfstats_outside_giab = 'data/NA24143.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}.outside_GIAB.GQ30.vcf.stats',
        NA24143_30x_runtime = 'data/NA24143.{aj_trio_build}/variants/longshot.pacbio.{aj_trio_aligner}.30x._/{chrom}.vcf.runtime',
        NA24143_30_cov = 'data/NA24143.{aj_trio_build}/aligned_reads/pacbio/pacbio.{aj_trio_aligner}.all.30x.bam.median_coverage',
    output:
        table = 'data/output/four_GIAB_genomes_table_extended.aj_trio_{aj_trio_build}_{aj_trio_aligner}.{chrom}.tex'
    run:
        NA12878_30x_table_files = genomes_table_files(input.NA12878_30x_vcfeval,
                                                  input.NA12878_30x_vcfstats_genome,
                                                  input.NA12878_30x_vcfstats_outside_giab,
                                                  input.NA12878_30x_runtime)

        NA12878_44x_table_files = genomes_table_files(input.NA12878_44x_vcfeval,
                                                  input.NA12878_44x_vcfstats_genome,
                                                  input.NA12878_44x_vcfstats_outside_giab,
                                                  input.NA12878_44x_runtime)

        NA24385_30x_table_files = genomes_table_files(input.NA24385_30x_vcfeval,
                                                  input.NA24385_30x_vcfstats_genome,
                                                  input.NA24385_30x_vcfstats_outside_giab,
                                                  input.NA24385_30x_runtime)

        NA24385_40x_table_files = genomes_table_files(input.NA24385_40x_vcfeval,
                                                  input.NA24385_40x_vcfstats_genome,
                                                  input.NA24385_40x_vcfstats_outside_giab,
                                                  input.NA24385_40x_runtime)

        NA24385_50x_table_files = genomes_table_files(input.NA24385_50x_vcfeval,
                                                  input.NA24385_50x_vcfstats_genome,
                                                  input.NA24385_50x_vcfstats_outside_giab,
                                                  input.NA24385_50x_runtime)

        NA24385_69x_table_files = genomes_table_files(input.NA24385_69x_vcfeval,
                                                  input.NA24385_69x_vcfstats_genome,
                                                  input.NA24385_69x_vcfstats_outside_giab,
                                                  input.NA24385_69x_runtime)

        NA24149_32x_table_files = genomes_table_files(input.NA24149_32x_vcfeval,
                                                  input.NA24149_32x_vcfstats_genome,
                                                  input.NA24149_32x_vcfstats_outside_giab,
                                                  input.NA24149_32x_runtime)

        NA24143_30x_table_files = genomes_table_files(input.NA24143_30x_vcfeval,
                                                  input.NA24143_30x_vcfstats_genome,
                                                  input.NA24143_30x_vcfstats_outside_giab,
                                                  input.NA24143_30x_runtime)
        covs = [parse_int_file(input.NA12878_30_cov),
         parse_int_file(input.NA12878_44_cov),
         parse_int_file(input.NA24385_30_cov),
         parse_int_file(input.NA24385_40_cov),
         parse_int_file(input.NA24385_50_cov),
         parse_int_file(input.NA24385_69_cov),
         parse_int_file(input.NA24149_32_cov),
         parse_int_file(input.NA24143_30_cov)]
        ptf.make_table_4_genomes_extended(NA12878_30x_table_files, NA12878_44x_table_files,
                                 NA24385_30x_table_files,
                                 NA24385_40x_table_files, NA24385_50x_table_files,
                                 NA24385_69x_table_files,
                                 NA24149_32x_table_files,
                                 NA24143_30x_table_files,
                                 covs,
                                 covs,
                                 output.table)

rule vcf_stats:
    params: job_name = lambda wildcards: 'vcf_stats.{}'.format(str(wildcards.x).replace("/", "."))
    input:
        vcfgz = '{x}.vcf.gz',
        tbi = '{x}.vcf.gz.tbi'
    output:
        stats = '{x}.vcf.stats'
    shell: '{RTGTOOLS} RTG_MEM=12g vcfstats {input.vcfgz} > {output.stats}'

rule filter_vcf_outside_GIAB:
    params: job_name = 'filter_vcf_outside_GIAB.{dataset}.{calls_name}.{chrom}.GQ{GQ}'
    input:  vcfgz = 'data/{dataset}/variants/{calls_name}/{chrom}.vcf.gz',
            region_filter = 'data/{dataset}/variants/ground_truth/region_filter.bed'
    output: vcfgz = 'data/{dataset}/variants/{calls_name}/{chrom}.outside_GIAB.GQ{GQ}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-bed={input.region_filter} -g {wildcards.GQ} -i {input.vcfgz} -o {output.vcfgz}'

rule filter_vcf_GQ:
    params: job_name = 'filter_vcf_GQ.{dataset}.{calls_name}.{chrom}.GQ{GQ}'
    input:  vcfgz = 'data/{dataset}/variants/{calls_name}/{chrom}.vcf.gz'
    output: vcfgz = 'data/{dataset}/variants/{calls_name}/{chrom}.GQ{GQ,\d+}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --keep-filter="PASS","." -g {wildcards.GQ} -i {input.vcfgz} -o {output.vcfgz}'

sim_covs = [20,30,40,60]
rule plot_simulation_pr_bars:
    params: job_name = 'plot_simulation_pr_bars.{chrom}.GQ{GQ}'
    input:
        longshot_genome = expand('data/simulation.1000g/vcfeval/longshot.pacbio.blasr.{c}x._/{{chrom}}',c=sim_covs),
        illumina_genome = expand('data/simulation.1000g/vcfeval/freebayes.illumina.aligned.{c}x.filtered/{{chrom}}',c=sim_covs),
        longshot_segdup = expand('data/simulation.1000g/vcfeval_segdup/longshot.pacbio.blasr.{c}x._/{{chrom}}',c=sim_covs),
        illumina_segdup = expand('data/simulation.1000g/vcfeval_segdup/freebayes.illumina.aligned.{c}x.filtered/{{chrom}}',c=sim_covs),
    output:
        png = 'data/plots/simulation_pr_barplot_genome_vs_segdup.{chrom}.GQ{GQ}.png'
    run:
        ptf.plot_precision_recall_bars_simulation(
            [x for x in input.longshot_genome],
            [x for x in input.illumina_genome],
            [x for x in input.longshot_segdup],
            [x for x in input.illumina_segdup],
            float(wildcards.GQ),
            [str(cov)+r'$\times$' for cov in sim_covs],
            output.png)

rule plot_simulation_pr_bars_extended:
    params: job_name = 'plot_simulation_pr_bars_extended.{chrom}.GQ{GQ}'
    input:
        longshot_ngmlr_genome = expand('data/simulation.1000g/vcfeval/longshot.pacbio.ngmlr.{c}x._/{{chrom}}',c=sim_covs),
        longshot_minimap2_genome = expand('data/simulation.1000g/vcfeval/longshot.pacbio.minimap2.{c}x._/{{chrom}}',c=sim_covs),
        longshot_bwamem_genome = expand('data/simulation.1000g/vcfeval/longshot.pacbio.bwa.{c}x._/{{chrom}}',c=sim_covs),
        longshot_blasr_genome = expand('data/simulation.1000g/vcfeval/longshot.pacbio.blasr.{c}x._/{{chrom}}',c=sim_covs),
        illumina_genome = expand('data/simulation.1000g/vcfeval/freebayes.illumina.aligned.{c}x.filtered/{{chrom}}',c=sim_covs),
        longshot_ngmlr_segdup = expand('data/simulation.1000g/vcfeval_segdup/longshot.pacbio.ngmlr.{c}x._/{{chrom}}',c=sim_covs),
        longshot_minimap2_segdup = expand('data/simulation.1000g/vcfeval_segdup/longshot.pacbio.minimap2.{c}x._/{{chrom}}',c=sim_covs),
        longshot_bwamem_segdup = expand('data/simulation.1000g/vcfeval_segdup/longshot.pacbio.bwa.{c}x._/{{chrom}}',c=sim_covs),
        longshot_blasr_segdup = expand('data/simulation.1000g/vcfeval_segdup/longshot.pacbio.blasr.{c}x._/{{chrom}}',c=sim_covs),
        illumina_segdup = expand('data/simulation.1000g/vcfeval_segdup/freebayes.illumina.aligned.{c}x.filtered/{{chrom}}',c=sim_covs),
    output:
        png = 'data/plots/simulation_pr_barplot_genome_vs_segdup_extended.{chrom}.GQ{GQ}.png'
    run:
        ptf.plot_precision_recall_bars_simulation_extended(
            pacbio_ngmlr_dirlist_genome=[x for x in input.longshot_ngmlr_genome],
            pacbio_minimap2_dirlist_genome=[x for x in input.longshot_minimap2_genome],
            pacbio_bwamem_dirlist_genome=[x for x in input.longshot_bwamem_genome],
            pacbio_blasr_dirlist_genome=[x for x in input.longshot_blasr_genome],
            illumina_dirlist_genome=[x for x in input.illumina_genome],
            pacbio_ngmlr_dirlist_segdup=[x for x in input.longshot_ngmlr_segdup],
            pacbio_minimap2_dirlist_segdup=[x for x in input.longshot_minimap2_segdup],
            pacbio_bwamem_dirlist_segdup=[x for x in input.longshot_bwamem_segdup],
            pacbio_blasr_dirlist_segdup=[x for x in input.longshot_blasr_segdup],
            illumina_dirlist_segdup=[x for x in input.illumina_segdup],
            gq_cutoff=float(wildcards.GQ),
            labels=[str(cov)+r'$\times$' for cov in sim_covs],
            output_file=output.png
        )

rule plot_precision_recall_bars_NA12878_AJ_Trio:
    params: job_name = 'plot_precision_recall_bars_NA12878_AJ_Trio.{AJ_trio_aligner}.{AJ_trio_build}',
    input: pacbio_dirlist_NA12878 = expand('data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.{cov}x._/all',cov=[30,44]),
           illumina_NA12878 = 'data/NA12878.1000g/vcfeval/freebayes.illumina.aligned.30x.filtered/all',
           pacbio_dirlist_NA24385 = expand('data/NA24385.{{AJ_trio_build}}/vcfeval/longshot.pacbio.{{AJ_trio_aligner}}.{cov}x._/all',cov=[30,40,50,69]),
           illumina_NA24385 = 'data/NA24385.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/all',
           pacbio_dirlist_NA24149 = expand('data/NA24149.{{AJ_trio_build}}/vcfeval/longshot.pacbio.{{AJ_trio_aligner}}.{cov}x._/all',cov=[32]),
           illumina_NA24149 = 'data/NA24149.hg38/vcfeval/freebayes.illumina.aligned.34x.filtered/all',
           pacbio_dirlist_NA24143 = expand('data/NA24143.{{AJ_trio_build}}/vcfeval/longshot.pacbio.{{AJ_trio_aligner}}.{cov}x._/all',cov=[30]),
           illumina_NA24143 = 'data/NA24143.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/all',
           NA12878_pb30_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
           NA12878_pb44_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.median_coverage',
           NA12878_il30_cov = 'data/NA12878.1000g/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           NA24385_pb30_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
           NA24385_pb40_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.40x.bam.median_coverage',
           NA24385_pb50_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.50x.bam.median_coverage',
           NA24385_pb69_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.69x.bam.median_coverage',
           NA24385_il30_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           NA24149_pb32_cov = 'data/NA24149.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.32x.bam.median_coverage',
           NA24149_il34_cov = 'data/NA24149.{AJ_trio_build}/aligned_reads/illumina/illumina.aligned.all.34x.bam.median_coverage',
           NA24143_pb30_cov = 'data/NA24143.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
           NA24143_il30_cov = 'data/NA24143.{AJ_trio_build}/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
    output: png = 'data/plots/fig3_precision_recall_bars_NA12878_AJ_Trio.{AJ_trio_aligner}.{AJ_trio_build}.png'
    run:
        ptf.plot_precision_recall_bars_NA12878_AJ_Trio(pacbio_dirlist_NA12878=list(input.pacbio_dirlist_NA12878),
                                                       illumina_dirlist_NA12878=[input.illumina_NA12878],
                                                       pacbio_dirlist_NA24385=list(input.pacbio_dirlist_NA24385),
                                                       illumina_dirlist_NA24385=[input.illumina_NA24385],
                                                       pacbio_dirlist_NA24149=list(input.pacbio_dirlist_NA24149),
                                                       illumina_dirlist_NA24149=[input.illumina_NA24149],
                                                       pacbio_dirlist_NA24143=list(input.pacbio_dirlist_NA24143),
                                                       illumina_dirlist_NA24143=[input.illumina_NA24143],
                                                       gq_cutoffs_NA12878=[parse_int_file(x) for x in
                                                                     [input.NA12878_pb30_cov, input.NA12878_pb44_cov]],
                                                       gq_cutoffs_NA24385=[parse_int_file(x) for x in
                                                                     [input.NA24385_pb30_cov, input.NA24385_pb40_cov,
                                                                      input.NA24385_pb50_cov, input.NA24385_pb69_cov]],
                                                       gq_cutoffs_NA24149=[parse_int_file(input.NA24149_pb32_cov)],
                                                       gq_cutoffs_NA24143=[parse_int_file(input.NA24143_pb30_cov)],
                                                       gq_cutoff_illumina=50,
                                                       labels = [str(parse_int_file(x))+r'$\times$' for x in
                                                                     [input.NA12878_pb30_cov, input.NA12878_pb44_cov,
                                                                      input.NA12878_il30_cov,
                                                                      input.NA24385_pb30_cov, input.NA24385_pb40_cov,
                                                                      input.NA24385_pb50_cov, input.NA24385_pb69_cov,
                                                                      input.NA24385_il30_cov,
                                                                      input.NA24149_pb32_cov, input.NA24149_il34_cov,
                                                                      input.NA24143_pb30_cov, input.NA24143_il30_cov,
                                                                      ]
                                                                ],
                                                       output_file=output.png)


rule plot_haplotyping_results:
    params: job_name = 'plot_haplotyping_results',
    input: longshot_NA12878_44x = 'data/NA12878.1000g/longshot_haplotypes/hap_statistics/longshot.pacbio.blasr.44x._.all.p',
           longshot_NA24385_69x = 'data/NA24385.hg38/longshot_haplotypes/hap_statistics/longshot.pacbio.blasr.69x._.all.p',
           hapcut2_NA12878_44x = 'data/NA12878.1000g/HapCUT2_haplotypes/hap_statistics/illumina.30x.pacbio.blasr.44x.all.p',
           hapcut2_NA24385_69x = 'data/NA24385.hg38/HapCUT2_haplotypes/hap_statistics/illumina.30x.pacbio.blasr.69x.all.p',
           NA12878_44_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.median_coverage',
           NA24385_69_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam.median_coverage',
    output: png = 'data/plots/haplotyping_results_barplot.png'
    run:
        ptf.plot_haplotyping_results(longshot_errs=[input.longshot_NA12878_44x,
                                                  input.longshot_NA24385_69x],
                                     hapcut2_errs=[input.hapcut2_NA12878_44x,
                                                  input.hapcut2_NA24385_69x],
                                     median_covs=[parse_int_file(input.NA12878_44_cov),
                                     parse_int_file(input.NA24385_69_cov)],
                                     output_file=output.png)

rule plot_precision_recall_bars_NA12878_AJ_Trio_with_haplotyping_results:
    params: job_name = 'plot_precision_recall_bars_NA12878_AJ_Trio_with_haplotyping_results.{AJ_trio_aligner}.{AJ_trio_build}',
    input: pacbio_dirlist_NA12878 = expand('data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.{cov}x._/all',cov=[30,44]),
           illumina_NA12878 = 'data/NA12878.1000g/vcfeval/freebayes.illumina.aligned.30x.filtered/all',
           pacbio_dirlist_NA24385 = expand('data/NA24385.{{AJ_trio_build}}/vcfeval/longshot.pacbio.{{AJ_trio_aligner}}.{cov}x._/all',cov=[30,40,50,69]),
           illumina_NA24385 = 'data/NA24385.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/all',
           pacbio_dirlist_NA24149 = expand('data/NA24149.{{AJ_trio_build}}/vcfeval/longshot.pacbio.{{AJ_trio_aligner}}.{cov}x._/all',cov=[32]),
           illumina_NA24149 = 'data/NA24149.hg38/vcfeval/freebayes.illumina.aligned.34x.filtered/all',
           pacbio_dirlist_NA24143 = expand('data/NA24143.{{AJ_trio_build}}/vcfeval/longshot.pacbio.{{AJ_trio_aligner}}.{cov}x._/all',cov=[30]),
           illumina_NA24143 = 'data/NA24143.hg38/vcfeval/freebayes.illumina.aligned.30x.filtered/all',
           NA12878_pb30_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
           NA12878_pb44_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.median_coverage',
           NA12878_il30_cov = 'data/NA12878.1000g/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           NA24385_pb30_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
           NA24385_pb40_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.40x.bam.median_coverage',
           NA24385_pb50_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.50x.bam.median_coverage',
           NA24385_pb69_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.69x.bam.median_coverage',
           NA24385_il30_cov = 'data/NA24385.{AJ_trio_build}/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           NA24149_pb32_cov = 'data/NA24149.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.32x.bam.median_coverage',
           NA24149_il34_cov = 'data/NA24149.{AJ_trio_build}/aligned_reads/illumina/illumina.aligned.all.34x.bam.median_coverage',
           NA24143_pb30_cov = 'data/NA24143.{AJ_trio_build}/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
           NA24143_il30_cov = 'data/NA24143.{AJ_trio_build}/aligned_reads/illumina/illumina.aligned.all.30x.bam.median_coverage',
           longshot_NA12878_30x_hap = 'data/NA12878.1000g/longshot_haplotypes/hap_statistics/longshot.pacbio.blasr.30x._.all.p',
           longshot_NA24385_30x_hap = 'data/NA24385.hg38/longshot_haplotypes/hap_statistics/longshot.pacbio.blasr.30x._.all.p',
           illumina_NA12878_30x_hap = 'data/NA12878.1000g/HapCUT2_haplotypes/hap_statistics/illumina.30x.illumina.aligned.30x.all.p',
           illumina_NA24385_30x_hap = 'data/NA24385.hg38/HapCUT2_haplotypes/hap_statistics/illumina.30x.illumina.aligned.30x.all.p'
    output: png = 'data/plots/fig3_precision_recall_bars_NA12878_AJ_Trio_with_haplotyping_results.{AJ_trio_aligner}.{AJ_trio_build}.png'
    run:
        ptf.plot_precision_recall_bars_NA12878_AJ_Trio_with_haplotyping_results(pacbio_dirlist_NA12878=list(input.pacbio_dirlist_NA12878),
                                                       illumina_dirlist_NA12878=[input.illumina_NA12878],
                                                       pacbio_dirlist_NA24385=list(input.pacbio_dirlist_NA24385),
                                                       illumina_dirlist_NA24385=[input.illumina_NA24385],
                                                       pacbio_dirlist_NA24149=list(input.pacbio_dirlist_NA24149),
                                                       illumina_dirlist_NA24149=[input.illumina_NA24149],
                                                       pacbio_dirlist_NA24143=list(input.pacbio_dirlist_NA24143),
                                                       illumina_dirlist_NA24143=[input.illumina_NA24143],
                                                       gq_cutoffs_NA12878=[parse_int_file(x) for x in
                                                                     [input.NA12878_pb30_cov, input.NA12878_pb44_cov]],
                                                       gq_cutoffs_NA24385=[parse_int_file(x) for x in
                                                                     [input.NA24385_pb30_cov, input.NA24385_pb40_cov,
                                                                      input.NA24385_pb50_cov, input.NA24385_pb69_cov]],
                                                       gq_cutoffs_NA24149=[parse_int_file(input.NA24149_pb32_cov)],
                                                       gq_cutoffs_NA24143=[parse_int_file(input.NA24143_pb30_cov)],
                                                       gq_cutoff_illumina=50,
                                                       labels_variants = [str(parse_int_file(x))+r'$\times$' for x in
                                                                     [input.NA12878_pb30_cov, input.NA12878_pb44_cov,
                                                                      input.NA12878_il30_cov,
                                                                      input.NA24385_pb30_cov, input.NA24385_pb40_cov,
                                                                      input.NA24385_pb50_cov, input.NA24385_pb69_cov,
                                                                      input.NA24385_il30_cov,
                                                                      input.NA24149_pb32_cov, input.NA24149_il34_cov,
                                                                      input.NA24143_pb30_cov, input.NA24143_il30_cov,
                                                                      ]
                                                                ],
                                                        labels_hap = [str(parse_int_file(x))+r'$\times$' for x in
                                                                      [input.NA12878_pb30_cov,
                                                                       input.NA12878_il30_cov,
                                                                       input.NA24385_pb30_cov,
                                                                       input.NA24385_il30_cov]],
                                                        longshot_errs=[input.longshot_NA12878_30x_hap,
                                                                       input.longshot_NA24385_30x_hap],
                                                        illumina_errs=[input.illumina_NA12878_30x_hap,
                                                                      input.illumina_NA24385_30x_hap],
                                                        output_file=output.png)

rule plot_fp_near_indel:
    params: job_name = 'plot_fp_near_indel',
    input: vcfevals = expand('data/NA24385.hg38/vcfeval/longshot.pacbio.blasr.{cov}x._/all',cov=[20,30,40,50,69]),
           ground_truth_vcfgz = 'data/NA24385.hg38/variants/ground_truth/ground_truth.vcf.gz',
           ground_truth_ix = 'data/NA24385.hg38/variants/ground_truth/ground_truth.vcf.gz.tbi',
           fixed_gq_VCFstats = expand('data/NA24385.hg38/variants/longshot.pacbio.blasr.{cov}x._/all.GQ30.PASS.SNPs_ONLY.vcf.stats',cov=[20,30,40,50,69]),
           scaled_gq_VCFstats = expand('data/NA24385.hg38/variants/longshot.pacbio.blasr.{cov}x._/all.GQ{cov}.PASS.SNPs_ONLY.vcf.stats',cov=[20,30,40,50,69]),
           NA24385_covs = expand('data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.{c}x.bam.median_coverage',c=[20,30,40,50,69]),
    output: png = 'data/plots/plot_fp_near_indel.NA24385.hg38.png'
    run:
        fp_vcfs = [os.path.join(f, 'fp.vcf.gz') for f in input.vcfevals]
        med_covs = [parse_int_file(f) for f in input.NA24385_covs]
        ptf.plot_fp_near_indel(fp_vcfs, input.fixed_gq_VCFstats, input.scaled_gq_VCFstats, input.ground_truth_vcfgz, output.png,
                             cov=med_covs,fixed_gq=30, scaled_gqs=med_covs)


rule plot_actual_vs_effective_coverage:
    params: job_name = 'plot_actual_vs_effective_coverage.NA12878',
    input: vcfgz = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/1.vcf.gz'
    output: png = 'data/plots/actual_vs_effective_coverage.chr1.NA12878.44x.png'
    run:
        ptf.actual_to_effective_read_coverage_plot(input.vcfgz, output.png)

map_covs = [20,30,40]
TOTAL_BASES_HG19 = 2867437753 - (151100560 + 25653566 + 3675142) # ungapped genome length (non-N) minus (chrX+chrY+un) # https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
rule plot_mappability_bars:
    params: job_name = 'plot_mappability_bars',
    input: pacbio_mf_0_0 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.0/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_0_5 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.5/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_0_75 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.75/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_0_9 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.9/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_1_0 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac1.0/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_0 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.0/illumina.aligned.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_5 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.5/illumina.aligned.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_75 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.75/illumina.aligned.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_9 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.9/illumina.aligned.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_1_0 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac1.0/illumina.aligned.60x__all.map_count.txt',cov=map_covs),
    output: png = 'data/plots/plot_mappability_bars.simulation.1000g.png'
    run:
        ptf.plot_mappability_bars(pacbio_mf_0_0 = input.pacbio_mf_0_0,
                                  pacbio_mf_0_5 = input.pacbio_mf_0_5,
                                  pacbio_mf_0_75 = input.pacbio_mf_0_75,
                                  pacbio_mf_0_9 = input.pacbio_mf_0_9,
                                  pacbio_mf_1_0 = input.pacbio_mf_1_0,
                                  illumina_mf_0_0 = input.illumina_mf_0_0,
                                  illumina_mf_0_5 = input.illumina_mf_0_5,
                                  illumina_mf_0_75 = input.illumina_mf_0_75,
                                  illumina_mf_0_9 = input.illumina_mf_0_9,
                                  illumina_mf_1_0 = input.illumina_mf_1_0,
                                  output_file=output.png,
                                  total_nonN_bases=TOTAL_BASES_HG19)

rule combine_map_counts:
    params: job_name = 'combine_map_counts.{individual}.{build}.{tech}.{info}.mapq{mapq}.mincov{cov}.mapfrac{mapfrac}'
    input:  expand('data/{{individual}}.{{build}}/aligned_reads/{{tech}}/map_counts/mapq{{mapq}}_mincov{{cov}}_mapfrac{{mapfrac}}/{{info}}__{chrom}.map_count.txt', chrom=chroms)
    output: 'data/{individual}.{build}/aligned_reads/{tech}/map_counts/mapq{mapq}_mincov{cov}_mapfrac{mapfrac}/{info}__all.map_count.txt'
    run:
        total = 0
        for infile in input:
            with open(infile,'r') as inf:
                total += int(inf.readline().strip())

        with open(output[0],'w') as outf:
            print(total, file=outf)

rule generate_map_counts:
    params: job_name = 'generate_map_counts_{tech}.{aligner}.{individual}.{build}.mapq{mapq}.mincov{cov}.mapfrac{mapfrac}.{chrom}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.60x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.60x.bam.bai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: 'data/{individual}.{build}/aligned_reads/{tech}/map_counts/mapq{mapq,\d+}_mincov{cov,\d+}_mapfrac{mapfrac}/{tech}.{aligner}.60x__{chrom,(\d+)}.map_count.txt'
    run:
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        ref_fa = input.hs37d5 if wildcards.build == '1000g' else input.hg38
        # NA12878 is an exception, actually has hg19 chrom names instead of 1000g chrom names
        if wildcards.individual == 'NA12878' and wildcards.build == '1000g' and wildcards.tech=='pacbio':
            w_chrom = 'chr' + w_chrom
        shell('''
        {MAP_COUNTER} --chrom {w_chrom} --bam {input.bam} --ref {ref_fa} --min_cov {wildcards.cov} \
        --min_mapq {wildcards.mapq} --map_frac {wildcards.mapfrac} > {output}
        ''')

rule generate_aj_trio_table:
    params: job_name = 'generate_aj_trio_table.{aligner}'
    input: pacbio_row = expand('data/aj_trio.{{aligner}}/trio_shared_variant_sites/pacbio/mendelian/all.{region}.report.txt', region=['whole_genome','confident','nonconfident','segdup95'])
    output: table = 'data/output/pacbio_mendelian_concordance_table.{aligner}.tex'
    run:
        ptf.mendelian_concordance_table(input.pacbio_row, output.table)

rule analyze_variants:
    params: job_name = 'analyze_variants.{individual}.{build}.{aligner}.{cov}x.{chrom}.GQ{GQ}',
    input:  vcfeval = 'data/{individual}.{build}/vcfeval/longshot.pacbio.{aligner}.{cov}x._/{chrom}',
            pacbio_calls_vcfgz = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{cov}x._/{chrom}.GQ{GQ}.vcf.gz',
            pacbio_calls_ix = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{cov}x._/{chrom}.GQ{GQ}.vcf.gz.tbi',
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
    output: tex = 'data/output/variant_analysis_fp_fn_{individual,NA\d+}.{build,(1000g|hg38)}.{aligner,(blasr|minimap2)}.{cov,\d+}x.GQ{GQ,\d+}.{chrom,(\d+|all)}.tex'
    run:
        analyze_variants(chrom_name = chr_prefix(wildcards.chrom, wildcards.build),
                         pacbio_calls_vcfgz = input.pacbio_calls_vcfgz,
                         fp_calls_vcfgz = os.path.join(input.vcfeval,'fp.vcf.gz'),
                         fn_calls_vcfgz =  os.path.join(input.vcfeval,'fn.vcf.gz'),
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
    input:  sizes_file = 'genome_tracks/{build}.chrom.sizes.txt'
    output: vcf = 'data/{individual}.{build,(1000g|hg38)}/random_positions/{chrom}.whole_chrom.vcf',
            vcfgz = 'data/{individual}.{build,(1000g|hg38)}/random_positions/{chrom}.whole_chrom.vcf.gz'
    run:
        chrom_name = chr_prefix(wildcards.chrom, wildcards.build)
        chrlen = get_chr_len(chrom_name, input.sizes_file)
        generate_random_calls(chrom_name, chrlen, 100000, output.vcf)
        shell('bgzip -c {output.vcf} > {output.vcfgz}')

rule copy_ground_truth:
    params: job_name = 'copy_ground_truth',
    input:  'data/NA12878.1000g/variants/ground_truth/ground_truth.vcf.gz'
    output: 'data/NA12878.1000g/variants/ground_truth/all.vcf.gz'
    shell: 'cp {input} {output}'

rule make_venn_diagram_variants_outside_GIAB:
    params: job_name = 'make_venn_diagram_variants_outside_GIAB',
    input:  longshot = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/all.outside_GIAB.GQ44.vcf.gz',
            giab = 'data/NA12878.1000g/variants/ground_truth/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz',
            plat = 'data/NA12878.1000g/variants/platinum_genomes/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz'
    output: png = 'data/plots/NA12878_variants_outside_GIAB_confident_venn_diagram.png'
    run:
        ptf.make_venn_diagram_variants_outside_GIAB(input.longshot, input.giab, input.plat, output.png)

NA12878_PG_CONFIDENT_REGIONS_URL = 'https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/ConfidentRegions.bed.gz'
rule download_pg_conf_NA12878:
    params: job_name = 'download_pg_conf_NA12878',
    input: names = 'genome_tracks/names.txt'
    output: bed = 'data/NA12878.1000g/variants/platinum_genomes/confident_regions.bed',
    shell:
        '''
        wget {NA12878_PG_CONFIDENT_REGIONS_URL} -O - | \
        gunzip -c | \
        python3 filter_bed_chroms.py --remove_chr | \
        {BEDTOOLS} sort -g {input.names} > {output.bed}
        '''

rule filter_SNVs_outside_confident_regions_inside_PG_conf:
    params: job_name = 'filter_SNVs_outside_GIAB_confident_inside_PG_conf.{individual}.{build}.{info}.GQ{GQ}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz.tbi',
            bed   = 'data/{individual}.{build}/variants/platinum_genomes/confident_regions.bed'
    output: vcfgz = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz'
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter --include-bed={input.bed} \
        -i {input.vcfgz} -o {output.vcfgz}
        '''

rule make_dual_venn_diagram_variants_outside_GIAB:
    params: job_name = 'make_venn_diagram_variants_outside_GIAB',
    input:  longshot = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/all.GQ44.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz', #'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/all.outside_GIAB.GQ44.vcf.gz',
            giab = 'data/NA12878.1000g/variants/ground_truth/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz',
            plat = 'data/NA12878.1000g/variants/platinum_genomes/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz',
            longshot_PGconf = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/all.GQ44.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
            giab_PGconf = 'data/NA12878.1000g/variants/ground_truth/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
            plat_PGconf = 'data/NA12878.1000g/variants/platinum_genomes/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
    output: png = 'data/plots/NA12878_variants_outside_GIAB_confident_dual_venn_diagram.png'
    run:
        ptf.make_dual_venn_diagram_variants_outside_GIAB(input.longshot, input.giab, input.plat,
                                                    input.longshot_PGconf, input.giab_PGconf, input.plat_PGconf,
                                                    output.png)

rule make_venn_diagram_variants_outside_GIAB_inside_PG:
    params: job_name = 'make_venn_diagram_variants_outside_GIAB',
    input: longshot = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
           giab = 'data/NA12878.1000g/variants/ground_truth/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
           plat = 'data/NA12878.1000g/variants/platinum_genomes/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
    output: png = 'data/plots/NA12878_variants_outside_GIAB_confident_inside_PG_confident_venn_diagram.png'
    run:
        ptf.make_venn_diagram_variants_outside_GIAB(input.longshot, input.giab, input.plat, output.png)

rule get_longshot_PG_venndiagram_vcfs:
    params: job_name = 'get_longshot_PG_venndiagram_vcfs',
    input: longshot = 'data/NA12878.1000g/variants/longshot.pacbio.blasr.44x._/all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
           giab = 'data/NA12878.1000g/variants/ground_truth/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
           plat = 'data/NA12878.1000g/variants/platinum_genomes/all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
    output: longshot_intersect_plat = 'data/NA12878.1000g/variants/misc/INTERSECT_PG_longshot.pacbio.blasr.44x._.all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
            longshot_minus_plat = 'data/NA12878.1000g/variants/misc/MINUS_longshot.pacbio.blasr.44x._.PG.all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz',
            plat_minus_longshot = 'data/NA12878.1000g/variants/misc/MINUS_PG_longshot.pacbio.blasr.44x._.all.GQ45.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.inside_PG_confident.vcf.gz'
    shell:
        '''
        {BEDTOOLS} intersect -header -a {input.longshot} -b {input.plat} | bgzip -c > {output.longshot_intersect_plat}
        {BEDTOOLS} subtract -header -a {input.longshot} -b {input.plat} | bgzip -c > {output.longshot_minus_plat}
        {BEDTOOLS} subtract -header -b {input.longshot} -a {input.plat} | bgzip -c > {output.plat_minus_longshot}
        '''

rule make_table_filtered_indels:
    params: job_name = 'make_table_filtered_indels',
    input:  NA12878_cov = 'data/NA12878.1000g/aligned_reads/pacbio/pacbio.blasr.all.44x.bam.median_coverage',
            NA24385_cov = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.all.69x.bam.median_coverage',
            NA24149_cov = 'data/NA24149.hg38/aligned_reads/pacbio/pacbio.blasr.all.32x.bam.median_coverage',
            NA24143_cov = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.all.30x.bam.median_coverage',
            NA12878_vcfeval_nofilter = 'data/NA12878.1000g/vcfeval/longshot.pacbio.blasr.44x._/all',
            NA12878_vcfeval_filter =   'data/NA12878.1000g/vcfeval/known_indels_removed.longshot.pacbio.blasr.44x._/all',
            NA24385_vcfeval_nofilter = 'data/NA24385.hg38/vcfeval/longshot.pacbio.blasr.69x._/all',
            NA24385_vcfeval_filter =   'data/NA24385.hg38/vcfeval/known_indels_removed.longshot.pacbio.blasr.69x._/all',
            NA24149_vcfeval_nofilter = 'data/NA24149.hg38/vcfeval/longshot.pacbio.blasr.32x._/all',
            NA24149_vcfeval_filter =   'data/NA24149.hg38/vcfeval/known_indels_removed.longshot.pacbio.blasr.32x._/all',
            NA24143_vcfeval_nofilter = 'data/NA24143.hg38/vcfeval/longshot.pacbio.blasr.30x._/all',
            NA24143_vcfeval_filter =   'data/NA24143.hg38/vcfeval/known_indels_removed.longshot.pacbio.blasr.30x._/all'
    output: tex = 'data/output/prec_recall_table_known_indels_filtered.tex'
    run:
        median_covs = [parse_int_file(input.NA12878_cov),
                       parse_int_file(input.NA24385_cov),
                       parse_int_file(input.NA24149_cov),
                       parse_int_file(input.NA24143_cov)]
        ptf.make_table_filtered_indels(median_covs,
                                 input.NA12878_vcfeval_nofilter, input.NA12878_vcfeval_filter,
                                 input.NA24385_vcfeval_nofilter, input.NA24385_vcfeval_filter,
                                 input.NA24149_vcfeval_nofilter, input.NA24149_vcfeval_filter,
                                 input.NA24143_vcfeval_nofilter, input.NA24143_vcfeval_filter,
                                 [44,69,32,30], output.tex)
