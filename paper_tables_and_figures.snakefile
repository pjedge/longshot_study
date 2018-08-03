from paper_tables_and_figures import genomes_table_files

rule plot_true_indel_FP_fixed_GQ:
    params: job_name = 'plot_true_indel_FP_fixed_GQ',
    input:
        expand('data/output/variant_analysis_fp_fn__NA24385.hg38__GQ50__reaper.pacbio.blasr.{cov}x.-z__1.tex',cov=[20,30,40,50,69])

rule plot_pr_curve_impact_of_haplotyping:
    params: job_name = 'plot_pr_curve_impact_of_haplotyping',
            title = 'Impact of Haplotype Information on PacBio Variant Calling with Reaper'
    input:
        NA12878_hap = 'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA24385_hap = 'data/NA24385.1000g/vcfeval/reaper.pacbio.bwamem.69x.-z/{chrom}.done',
        NA24149_hap = 'data/NA24149.1000g/vcfeval/reaper.pacbio.bwamem.32x.-z/{chrom}.done',
        NA24143_hap = 'data/NA24143.1000g/vcfeval/reaper.pacbio.bwamem.30x.-z/{chrom}.done',
        NA12878_no_hap = 'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.44x.-n_-z/{chrom}.done',
        NA24385_no_hap = 'data/NA24385.1000g/vcfeval/reaper.pacbio.bwamem.69x.-n_-z/{chrom}.done',
        NA24149_no_hap = 'data/NA24149.1000g/vcfeval/reaper.pacbio.bwamem.32x.-n_-z/{chrom}.done',
        NA24143_no_hap = 'data/NA24143.1000g/vcfeval/reaper.pacbio.bwamem.30x.-n_-z/{chrom}.done'
    output:
        png = 'data/plots/effect_of_haplotyping.giab_individuals.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.NA12878_hap[:-5], input.NA12878_no_hap[:-5],
                          input.NA24385_hap[:-5], input.NA24385_no_hap[:-5],
                          input.NA24149_hap[:-5], input.NA24149_no_hap[:-5],
                          input.NA24143_hap[:-5], input.NA24143_no_hap[:-5]],
                         ['NA12878, 44x (hap)', 'NA12878, 44x (no hap)',
                          'NA24385, 69x (hap)', 'NA24385, 69x (no hap)',
                          'NA24149, 32x (hap)', 'NA24149, 32x (no hap)',
                          'NA24143, 30x (hap)', 'NA24143, 30x (no hap)'],
                           output.png,params.title,
                           colors=['#ff0000','#ff4f4f', # red
                                   '#02fc0a','#4af950', # green
                                   '#0400fc','#514efc', # blue
                                   '#ffff02','#ffff6b'], # yellow
                           xlim=(0.5,1.0),ylim=(0.95,1.0))

rule plot_pr_curve_NA24385_impact_of_haplotyping:
    params: job_name = 'plot_pr_curve_impact_of_haplotyping',
            title = 'Impact of Haplotype Information on PacBio Variant Calling with Reaper'
    input:
        NA24385_hap_20x = 'data/NA24385.1000g/vcfeval/reaper.pacbio.bwamem.20x.-z/{chrom}.done',
        NA24385_hap_69x = 'data/NA24385.1000g/vcfeval/reaper.pacbio.bwamem.69x.-z/{chrom}.done',
        NA24385_no_hap_20x = 'data/NA24385.1000g/vcfeval/reaper.pacbio.bwamem.20x.-n_-z/{chrom}.done',
        NA24385_no_hap_69x = 'data/NA24385.1000g/vcfeval/reaper.pacbio.bwamem.69x.-n_-z/{chrom}.done',
    output:
        png = 'data/plots/effect_of_haplotyping.NA24385.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.NA24385_hap_20x[:-5], input.NA24385_no_hap_20x[:-5],
                          input.NA24385_hap_69x[:-5], input.NA24385_no_hap_69x[:-5]],
                         ['20x SMS\n(haplotype assembly)', '20x SMS\n(no haplotype assembly)',
                          '69x SMS\n(haplotype assembly)', '69x SMS\n(no haplotype assembly)'],
                           output.png,None,
                           colors=['#ff0000','#ff4f4f', # red
                                   '#0400fc','#514efc'], # blue
                           xlim=(0.70,1.0),ylim=(0.975,1.0))

rule make_four_genomes_table:
    params: job_name = 'make_four_genomes_table.aj_trio_{aj_trio_build}_{aj_trio_aligner}.{chrom}'
    input:
        NA12878_vcfeval = 'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA12878_vcfstats_genome = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.-z/{chrom}.GQ44.vcf.stats',
        NA12878_vcfstats_outside_giab = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.-z/{chrom}.outside_GIAB.GQ44.vcf.stats',
        NA12878_runtime = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.-z/{chrom}.vcf.runtime',
        NA24385_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.done',
        NA24385_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.GQ69.vcf.stats',
        NA24385_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.outside_GIAB.GQ69.vcf.stats',
        NA24385_runtime = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.vcf.runtime',
        NA24149_vcfeval = 'data/NA24149.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.done',
        NA24149_vcfstats_genome = 'data/NA24149.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.GQ32.vcf.stats',
        NA24149_vcfstats_outside_giab = 'data/NA24149.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.outside_GIAB.GQ32.vcf.stats',
        NA24149_runtime = 'data/NA24149.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.vcf.runtime',
        NA24143_vcfeval = 'data/NA24143.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.done',
        NA24143_vcfstats_genome = 'data/NA24143.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.GQ30.vcf.stats',
        NA24143_vcfstats_outside_giab = 'data/NA24143.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.outside_GIAB.GQ30.vcf.stats',
        NA24143_runtime = 'data/NA24143.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.vcf.runtime',
    output:
        table = 'data/output/four_GIAB_genomes_table.aj_trio_{aj_trio_build}_{aj_trio_aligner}.{chrom}.tex'
    run:
        NA12878_table_files = genomes_table_files(input.NA12878_vcfeval[:-5],
                                                  input.NA12878_vcfstats_genome,
                                                  input.NA12878_vcfstats_outside_giab,
                                                  input.NA12878_runtime)

        NA24385_table_files = genomes_table_files(input.NA24385_vcfeval[:-5],
                                                  input.NA24385_vcfstats_genome,
                                                  input.NA24385_vcfstats_outside_giab,
                                                  input.NA24385_runtime)

        NA24149_table_files = genomes_table_files(input.NA24149_vcfeval[:-5],
                                                  input.NA24149_vcfstats_genome,
                                                  input.NA24149_vcfstats_outside_giab,
                                                  input.NA24149_runtime)

        NA24143_table_files = genomes_table_files(input.NA24143_vcfeval[:-5],
                                                  input.NA24143_vcfstats_genome,
                                                  input.NA24143_vcfstats_outside_giab,
                                                  input.NA24143_runtime)
        ptf.make_table_4_genomes(NA12878_table_files, NA24385_table_files,
                                 NA24149_table_files, NA24143_table_files,
                                 [44,69,32,30], output.table)

rule make_four_genomes_table_extended:
    params: job_name = 'make_four_genomes_table_extended.aj_trio_{aj_trio_build}_{aj_trio_aligner}.{chrom}'
    input:
        NA12878_30x_vcfeval = 'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.30x.-z/{chrom}.done',
        NA12878_30x_vcfstats_genome = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.30x.-z/{chrom}.GQ30.vcf.stats',
        NA12878_30x_vcfstats_outside_giab = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.30x.-z/{chrom}.outside_GIAB.GQ30.vcf.stats',
        NA12878_30x_runtime = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.30x.-z/{chrom}.vcf.runtime',
        NA12878_44x_vcfeval = 'data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA12878_44x_vcfstats_genome = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.-z/{chrom}.GQ44.vcf.stats',
        NA12878_44x_vcfstats_outside_giab = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.-z/{chrom}.outside_GIAB.GQ44.vcf.stats',
        NA12878_44x_runtime = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.-z/{chrom}.vcf.runtime',
        NA24385_20x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.20x.-z/{chrom}.done',
        NA24385_20x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.20x.-z/{chrom}.GQ20.vcf.stats',
        NA24385_20x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.20x.-z/{chrom}.outside_GIAB.GQ20.vcf.stats',
        NA24385_20x_runtime = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.20x.-z/{chrom}.vcf.runtime',
        NA24385_30x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.done',
        NA24385_30x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.GQ30.vcf.stats',
        NA24385_30x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.outside_GIAB.GQ30.vcf.stats',
        NA24385_30x_runtime = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.vcf.runtime',
        NA24385_40x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.40x.-z/{chrom}.done',
        NA24385_40x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.40x.-z/{chrom}.GQ40.vcf.stats',
        NA24385_40x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.40x.-z/{chrom}.outside_GIAB.GQ40.vcf.stats',
        NA24385_40x_runtime = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.40x.-z/{chrom}.vcf.runtime',
        NA24385_50x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.50x.-z/{chrom}.done',
        NA24385_50x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.50x.-z/{chrom}.GQ50.vcf.stats',
        NA24385_50x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.50x.-z/{chrom}.outside_GIAB.GQ50.vcf.stats',
        NA24385_50x_runtime = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.50x.-z/{chrom}.vcf.runtime',
        NA24385_69x_vcfeval = 'data/NA24385.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.done',
        NA24385_69x_vcfstats_genome = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.GQ69.vcf.stats',
        NA24385_69x_vcfstats_outside_giab = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.outside_GIAB.GQ69.vcf.stats',
        NA24385_69x_runtime = 'data/NA24385.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.69x.-z/{chrom}.vcf.runtime',
        NA24149_32x_vcfeval = 'data/NA24149.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.done',
        NA24149_32x_vcfstats_genome = 'data/NA24149.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.GQ32.vcf.stats',
        NA24149_32x_vcfstats_outside_giab = 'data/NA24149.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.outside_GIAB.GQ32.vcf.stats',
        NA24149_32x_runtime = 'data/NA24149.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.32x.-z/{chrom}.vcf.runtime',
        NA24143_30x_vcfeval = 'data/NA24143.{aj_trio_build}/vcfeval/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.done',
        NA24143_30x_vcfstats_genome = 'data/NA24143.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.GQ30.vcf.stats',
        NA24143_30x_vcfstats_outside_giab = 'data/NA24143.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.outside_GIAB.GQ30.vcf.stats',
        NA24143_30x_runtime = 'data/NA24143.{aj_trio_build}/variants/reaper.pacbio.{aj_trio_aligner}.30x.-z/{chrom}.vcf.runtime',
    output:
        table = 'data/output/four_GIAB_genomes_table_extended.aj_trio_{aj_trio_build}_{aj_trio_aligner}.{chrom}.tex'
    run:
        NA12878_30x_table_files = genomes_table_files(input.NA12878_30x_vcfeval[:-5],
                                                  input.NA12878_30x_vcfstats_genome,
                                                  input.NA12878_30x_vcfstats_outside_giab,
                                                  input.NA12878_30x_runtime)

        NA12878_44x_table_files = genomes_table_files(input.NA12878_44x_vcfeval[:-5],
                                                  input.NA12878_44x_vcfstats_genome,
                                                  input.NA12878_44x_vcfstats_outside_giab,
                                                  input.NA12878_44x_runtime)

        NA24385_20x_table_files = genomes_table_files(input.NA24385_20x_vcfeval[:-5],
                                                  input.NA24385_20x_vcfstats_genome,
                                                  input.NA24385_20x_vcfstats_outside_giab,
                                                  input.NA24385_20x_runtime)

        NA24385_30x_table_files = genomes_table_files(input.NA24385_30x_vcfeval[:-5],
                                                  input.NA24385_30x_vcfstats_genome,
                                                  input.NA24385_30x_vcfstats_outside_giab,
                                                  input.NA24385_30x_runtime)

        NA24385_40x_table_files = genomes_table_files(input.NA24385_40x_vcfeval[:-5],
                                                  input.NA24385_40x_vcfstats_genome,
                                                  input.NA24385_40x_vcfstats_outside_giab,
                                                  input.NA24385_40x_runtime)

        NA24385_50x_table_files = genomes_table_files(input.NA24385_50x_vcfeval[:-5],
                                                  input.NA24385_50x_vcfstats_genome,
                                                  input.NA24385_50x_vcfstats_outside_giab,
                                                  input.NA24385_50x_runtime)

        NA24385_69x_table_files = genomes_table_files(input.NA24385_69x_vcfeval[:-5],
                                                  input.NA24385_69x_vcfstats_genome,
                                                  input.NA24385_69x_vcfstats_outside_giab,
                                                  input.NA24385_69x_runtime)

        NA24149_32x_table_files = genomes_table_files(input.NA24149_32x_vcfeval[:-5],
                                                  input.NA24149_32x_vcfstats_genome,
                                                  input.NA24149_32x_vcfstats_outside_giab,
                                                  input.NA24149_32x_runtime)

        NA24143_30x_table_files = genomes_table_files(input.NA24143_30x_vcfeval[:-5],
                                                  input.NA24143_30x_vcfstats_genome,
                                                  input.NA24143_30x_vcfstats_outside_giab,
                                                  input.NA24143_30x_runtime)
        ptf.make_table_4_genomes_extended(NA12878_30x_table_files, NA12878_44x_table_files,
                                 NA24385_20x_table_files, NA24385_30x_table_files,
                                 NA24385_40x_table_files, NA24385_50x_table_files,
                                 NA24385_69x_table_files,
                                 NA24149_32x_table_files,
                                 NA24143_30x_table_files,
                                 [30,44,20,30,40,50,69,32,30], output.table)

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
        reaper_genome = expand('data/simulation.1000g/vcfeval/reaper.pacbio.blasr.{c}x.-z/{{chrom}}.done',c=sim_covs),
        illumina_genome = expand('data/simulation.1000g/vcfeval/illumina_{c}x.filtered/{{chrom}}.done',c=sim_covs),
        reaper_segdup = expand('data/simulation.1000g/vcfeval_segdup/reaper.pacbio.blasr.{c}x.-z/{{chrom}}.done',c=sim_covs),
        illumina_segdup = expand('data/simulation.1000g/vcfeval_segdup/illumina_{c}x.filtered/{{chrom}}.done',c=sim_covs),
    output:
        png = 'data/plots/simulation_pr_barplot_genome_vs_segdup.{chrom}.GQ{GQ}.png'
    run:
        ptf.plot_precision_recall_bars_simulation(
            [x[:-5] for x in input.reaper_genome],
            [x[:-5] for x in input.illumina_genome],
            [x[:-5] for x in input.reaper_segdup],
            [x[:-5] for x in input.illumina_segdup],
            float(wildcards.GQ),
            sim_covs,
            output.png)

rule plot_simulation_pr_bars_extended:
    params: job_name = 'plot_simulation_pr_bars_extended.{chrom}.GQ{GQ}'
    input:
        reaper_ngmlr_genome = expand('data/simulation.1000g/vcfeval/reaper.pacbio.ngmlr.{c}x.-z/{{chrom}}.done',c=sim_covs),
        reaper_minimap2_genome = expand('data/simulation.1000g/vcfeval/reaper.pacbio.minimap2.{c}x.-z/{{chrom}}.done',c=sim_covs),
        reaper_bwamem_genome = expand('data/simulation.1000g/vcfeval/reaper.pacbio.bwa.{c}x.-z/{{chrom}}.done',c=sim_covs),
        reaper_blasr_genome = expand('data/simulation.1000g/vcfeval/reaper.pacbio.blasr.{c}x.-z/{{chrom}}.done',c=sim_covs),
        illumina_genome = expand('data/simulation.1000g/vcfeval/illumina_{c}x.filtered/{{chrom}}.done',c=sim_covs),
        reaper_ngmlr_segdup = expand('data/simulation.1000g/vcfeval_segdup/reaper.pacbio.ngmlr.{c}x.-z/{{chrom}}.done',c=sim_covs),
        reaper_minimap2_segdup = expand('data/simulation.1000g/vcfeval_segdup/reaper.pacbio.minimap2.{c}x.-z/{{chrom}}.done',c=sim_covs),
        reaper_bwamem_segdup = expand('data/simulation.1000g/vcfeval_segdup/reaper.pacbio.bwa.{c}x.-z/{{chrom}}.done',c=sim_covs),
        reaper_blasr_segdup = expand('data/simulation.1000g/vcfeval_segdup/reaper.pacbio.blasr.{c}x.-z/{{chrom}}.done',c=sim_covs),
        illumina_segdup = expand('data/simulation.1000g/vcfeval_segdup/illumina_{c}x.filtered/{{chrom}}.done',c=sim_covs),
    output:
        png = 'data/plots/simulation_pr_barplot_genome_vs_segdup_extended.{chrom}.GQ{GQ}.png'
    run:
        ptf.plot_precision_recall_bars_simulation_extended(
            pacbio_ngmlr_dirlist_genome=[x[:-5] for x in input.reaper_ngmlr_genome],
            pacbio_minimap2_dirlist_genome=[x[:-5] for x in input.reaper_minimap2_genome],
            pacbio_bwamem_dirlist_genome=[x[:-5] for x in input.reaper_bwamem_genome],
            pacbio_blasr_dirlist_genome=[x[:-5] for x in input.reaper_blasr_genome],
            illumina_dirlist_genome=[x[:-5] for x in input.illumina_genome],
            pacbio_ngmlr_dirlist_segdup=[x[:-5] for x in input.reaper_ngmlr_segdup],
            pacbio_minimap2_dirlist_segdup=[x[:-5] for x in input.reaper_minimap2_segdup],
            pacbio_bwamem_dirlist_segdup=[x[:-5] for x in input.reaper_bwamem_segdup],
            pacbio_blasr_dirlist_segdup=[x[:-5] for x in input.reaper_blasr_segdup],
            illumina_dirlist_segdup=[x[:-5] for x in input.illumina_segdup],
            gq_cutoff=float(wildcards.GQ),
            labels=sim_covs,
            output_file=output.png
        )

rule plot_depth_of_mapped_vs_breadth:
    params: job_name = 'plot_depth_of_mapped_vs_breadth.simulation.30x',
    input: illumina_mapq0 = 'data/simulation.1000g/aligned_reads/illumina/genomecov_histograms_mapq0/illumina.30x__all.txt',
           illumina_mapq10 = 'data/simulation.1000g/aligned_reads/illumina/genomecov_histograms_mapq10/illumina.30x__all.txt',
           illumina_mapq20 = 'data/simulation.1000g/aligned_reads/illumina/genomecov_histograms_mapq20/illumina.30x__all.txt',
           illumina_mapq30 = 'data/simulation.1000g/aligned_reads/illumina/genomecov_histograms_mapq30/illumina.30x__all.txt',
           pacbio_mapq0 = 'data/simulation.1000g/aligned_reads/pacbio/genomecov_histograms_mapq0/pacbio.blasr.all.30x__all.txt',
           pacbio_mapq10 = 'data/simulation.1000g/aligned_reads/pacbio/genomecov_histograms_mapq10/pacbio.blasr.all.30x__all.txt',
           pacbio_mapq20 = 'data/simulation.1000g/aligned_reads/pacbio/genomecov_histograms_mapq20/pacbio.blasr.all.30x__all.txt',
           pacbio_mapq30 = 'data/simulation.1000g/aligned_reads/pacbio/genomecov_histograms_mapq30/pacbio.blasr.all.30x__all.txt',
    output: png = 'data/plots/depth_vs_breadth_mappability.simulation.30x.png'
    run:

        inputs = [input.illumina_mapq0, input.illumina_mapq30, input.pacbio_mapq0, input.pacbio_mapq30]
        labels = ['Illumina, Mapq >= 0','Illumina, Mapq >= 30',
                  'PacBio, Mapq >= 0','PacBio, Mapq >= 30']

        #colors = ['#ff9999','#ff7070','#ff4242','#ff0000',
        #          '#9a99ff','#7775ff','#4744ff','#0400ff']

        colors = ['#ff9999', '#ff0000',
                  '#9a99ff', '#0400ff']

        #ptf.plot_depth_of_mapped_vs_breadth(list(input), labels, colors, output.png)
        ptf.plot_depth_of_mapped_vs_breadth(inputs, labels, colors, output.png)

rule plot_precision_recall_bars_NA12878_NA24385:
    params: job_name = 'plot_precision_recall_bars_NA12878_NA24385.{NA24385_aligner}.{NA24385_build}',
    input: pacbio_dirlist_NA24385 = expand('data/NA24385.{{NA24385_build}}/vcfeval/reaper.pacbio.{{NA24385_aligner}}.{cov}x.-z/all.done',cov=[20,30,40,50,69]),
           illumina_NA24385 = 'data/NA24385.1000g/vcfeval/illumina_30x.filtered/all.done',
           pacbio_dirlist_NA12878 = expand('data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.{cov}x.-z/all.done',cov=[30,44]),
           illumina_NA12878 = 'data/NA24385.1000g/vcfeval/illumina_30x.filtered/all.done',
    output: png = 'data/plots/fig3_precision_recall_bars_NA24385_NA12878.{NA24385_aligner}.{NA24385_build}.png'
    run:
        ptf.plot_precision_recall_bars_NA12878_NA24385(pacbio_dirlist_NA24385=[x[:-5] for x in input.pacbio_dirlist_NA24385],
                                                       illumina_dirlist_NA24385=[input.illumina_NA24385[:-5]],
                                                       pacbio_dirlist_NA12878=[x[:-5] for x in input.pacbio_dirlist_NA12878],
                                                       illumina_dirlist_NA12878=[input.illumina_NA12878[:-5]],
                                                       gq_cutoffs_NA24385=[20,30,40,50,69],gq_cutoffs_NA12878=[30,44],
                                                       gq_cutoff_illumina=50, output_file=output.png)


rule plot_haplotyping_results:
    params: job_name = 'plot_haplotyping_results',
    input: reaper_NA12878_30x = 'data/NA12878.1000g/reaper_haplotypes/hap_statistics/reaper.pacbio.blasr.30x.-z.all.p',
           reaper_NA12878_44x = 'data/NA12878.1000g/reaper_haplotypes/hap_statistics/reaper.pacbio.blasr.44x.-z.all.p',
           reaper_NA24385_69x = 'data/NA24385.hg38/reaper_haplotypes/hap_statistics/reaper.pacbio.blasr.69x.-z.all.p',
           reaper_NA24149_32x = 'data/NA24149.hg38/reaper_haplotypes/hap_statistics/reaper.pacbio.blasr.32x.-z.all.p',
           reaper_NA24143_30x = 'data/NA24143.hg38/reaper_haplotypes/hap_statistics/reaper.pacbio.blasr.30x.-z.all.p',
           hapcut2_NA12878_30x = 'data/NA12878.1000g/HapCUT2_haplotypes/hap_statistics/illumina.30x.pacbio.blasr.30x.all.p',
           hapcut2_NA12878_44x = 'data/NA12878.1000g/HapCUT2_haplotypes/hap_statistics/illumina.30x.pacbio.blasr.44x.all.p',
           hapcut2_NA24385_69x = 'data/NA24385.hg38/HapCUT2_haplotypes/hap_statistics/illumina.30x.pacbio.blasr.69x.all.p',
           hapcut2_NA24149_32x = 'data/NA24149.hg38/HapCUT2_haplotypes/hap_statistics/illumina.30x.pacbio.blasr.32x.all.p',
           hapcut2_NA24143_30x = 'data/NA24143.hg38/HapCUT2_haplotypes/hap_statistics/illumina.30x.pacbio.blasr.30x.all.p',
    output: png = 'data/plots/haplotyping_results_barplot.png'
    run:
        ptf.plot_haplotyping_results(reaper_errs=[input.reaper_NA12878_30x,
                                                  input.reaper_NA12878_44x,
                                                  input.reaper_NA24385_69x,
                                                  input.reaper_NA24149_32x,
                                                  input.reaper_NA24143_30x],
                                     hapcut2_errs=[input.hapcut2_NA12878_30x,
                                                  input.hapcut2_NA12878_44x,
                                                  input.hapcut2_NA24385_69x,
                                                  input.hapcut2_NA24149_32x,
                                                  input.hapcut2_NA24143_30x],
                                     output_file=output.png)


rule plot_precision_recall_bars_NA12878_NA24385_with_haplotypefree:
    params: job_name = 'plot_precision_recall_bars_NA12878_NA24385_with_haplotypefree.{NA24385_aligner}.{NA24385_build}',
    input: hap_dirlist_NA24385 = expand('data/NA24385.{{NA24385_build}}/vcfeval/reaper.pacbio.{{NA24385_aligner}}.{cov}x.-z/all.done',cov=[20,30,40,50,69]),
           nohap_dirlist_NA24385 = expand('data/NA24385.{{NA24385_build}}/vcfeval/reaper.pacbio.{{NA24385_aligner}}.{cov}x.-n_-z/all.done',cov=[20,30,40,50,69]),
           hap_dirlist_NA12878 = expand('data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.{cov}x.-z/all.done',cov=[30,44]),
           nohap_dirlist_NA12878 = expand('data/NA12878.1000g/vcfeval/reaper.pacbio.blasr.{cov}x.-n_-z/all.done',cov=[30,44]),
    output: png = 'data/plots/supp_fig3_with_haplotypefree_precision_recall_bars_NA24385_NA12878.{NA24385_aligner}.{NA24385_build}.png'
    run:
        ptf.plot_precision_recall_bars_NA12878_NA24385_with_haplotypefree(hap_dirlist_NA24385=[x[:-5] for x in input.hap_dirlist_NA24385],
                                                       nohap_dirlist_NA24385=[x[:-5] for x in input.nohap_dirlist_NA24385],
                                                       hap_dirlist_NA12878=[x[:-5] for x in input.hap_dirlist_NA12878],
                                                       nohap_dirlist_NA12878=[x[:-5] for x in input.nohap_dirlist_NA12878],
                                                       gq_cutoffs_NA24385=[20,30,40,50,69],gq_cutoffs_NA12878=[30,44],
                                                       output_file=output.png)

rule plot_fp_near_indel:
    params: job_name = 'plot_fp_near_indel',
    input: vcfevals = expand('data/NA24385.hg38/vcfeval/reaper.pacbio.blasr.{cov}x.-z/all.done',cov=[20,30,40,50,69]),
           ground_truth_vcfgz = 'data/NA24385.hg38/variants/ground_truth/ground_truth.vcf.gz',
           ground_truth_ix = 'data/NA24385.hg38/variants/ground_truth/ground_truth.vcf.gz.tbi',
           fixed_gq_VCFstats = expand('data/NA24385.hg38/variants/reaper.pacbio.blasr.{cov}x.-z/all.GQ30.PASS.SNPs_ONLY.vcf.stats',cov=[20,30,40,50,69]),
           scaled_gq_VCFstats = expand('data/NA24385.hg38/variants/reaper.pacbio.blasr.{cov}x.-z/all.GQ{cov}.PASS.SNPs_ONLY.vcf.stats',cov=[20,30,40,50,69]),
    output: png = 'data/plots/plot_fp_near_indel.NA24385.hg38.png'
    run:
        fp_vcfs = [os.path.join(f[:-5], 'fp.vcf.gz') for f in input.vcfevals]
        ptf.plot_fp_near_indel(fp_vcfs, input.fixed_gq_VCFstats, input.scaled_gq_VCFstats, input.ground_truth_vcfgz, output.png,
                             cov=[20,30,40,50,69],fixed_gq=30, scaled_gqs=[20,30,40,50,69])


rule plot_actual_vs_effective_coverage:
    params: job_name = 'plot_actual_vs_effective_coverage.NA12878',
    input: vcfgz = 'data/NA12878.1000g/variants/reaper.pacbio.blasr.44x.-z/1.vcf.gz'
    output: png = 'data/plots/actual_vs_effective_coverage.chr1.NA12878.44x.png'
    run:
        ptf.actual_to_effective_read_coverage_scatterplot(input.vcfgz, output.png)

map_covs = [10,20,30,40,50]
TOTAL_BASES_HG19 = 2867437753 - (151100560 + 25653566 + 3675142) # ungapped genome length (non-N) minus (chrX+chrY+un) # https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
rule plot_mappability_bars:
    params: job_name = 'plot_mappability_bars',
    input: pacbio_mf_0_0 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.0/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_0_5 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.5/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_0_75 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.75/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_0_9 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac0.9/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           pacbio_mf_1_0 = expand('data/simulation.1000g/aligned_reads/pacbio/map_counts/mapq30_mincov{cov}_mapfrac1.0/pacbio.blasr.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_0 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.0/illumina.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_5 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.5/illumina.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_75 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.75/illumina.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_0_9 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac0.9/illumina.60x__all.map_count.txt',cov=map_covs),
           illumina_mf_1_0 = expand('data/simulation.1000g/aligned_reads/illumina/map_counts/mapq30_mincov{cov}_mapfrac1.0/illumina.60x__all.map_count.txt',cov=map_covs),
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

rule generate_map_counts_illumina:
    params: job_name = 'generate_map_counts_illumina.{individual}.{build}.MAPQ{mapq}.mincov{cov}.mapfrac{mapfrac}.{chrom}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/illumina/illumina.60x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/illumina/illumina.60x.bam.bai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: 'data/{individual}.{build}/aligned_reads/illumina/map_counts/mapq{mapq,\d+}_mincov{cov,\d+}_mapfrac{mapfrac}/illumina.60x__{chrom,(\d+)}.map_count.txt'
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

rule generate_map_counts_pacbio:
    params: job_name = 'generate_map_counts_pacbio.{individual}.{build}.mapq{mapq}.mincov{cov}.mapfrac{mapfrac}.{chrom}'
    input:  bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.blasr.{chrom}.60x.bam',
            bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.blasr.{chrom}.60x.bam.bai',
            hs37d5    = 'data/genomes/hs37d5.fa',
            hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
            hg38 = 'data/genomes/hg38.fa',
            hg38_ix = 'data/genomes/hg38.fa.fai'
    output: 'data/{individual}.{build}/aligned_reads/pacbio/map_counts/mapq{mapq,\d+}_mincov{cov,\d+}_mapfrac{mapfrac}/pacbio.blasr.60x__{chrom,(\d+)}.map_count.txt'
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
