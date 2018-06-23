from paper_tables_and_figures import genomes_table_files

rule plot_pr_curve_NA12878_impact_of_haplotyping:
    params: job_name = 'plot_pr_curve_NA12878_impact_of_haplotyping',
            title = 'NA12878: Impact of Haplotype Information on PacBio Variant Calling with Reaper'
    input:
        NA12878_3_0 = 'data/NA12878/vcfeval_3.0/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA12878_no_hap = 'data/NA12878/vcfeval_no_haps/reaper.pacbio.blasr.44x.-z/{chrom}.done'
    output:
        png = 'data/plots/effect_of_haplotyping.NA12878.prec_recall_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.NA12878_3_0[:-5], input.NA12878_no_hap[:-5]],
                         ['NA12878, 44x (one round hap)',  'NA12878, 44x (no hap)'],
                           output.png,params.title,
                           colors=['b','y'],
                           xlim=(0.9,1.0),ylim=(0.99,1.0))

rule plot_pr_curve_varyQ:
    params: job_name = 'plot_pr_curve_varyQ',
            title = 'Impact of Haplotype Q Parameter'
    input:
        NA12878_30x_Q50 = 'data/NA12878/vcfeval/reaper.pacbio.blasr.30x.-z_-Q_50/{chrom}.done',
        NA12878_30x_Q20 = 'data/NA12878/vcfeval/reaper.pacbio.blasr.30x.-z_-Q_20/{chrom}.done',
        NA12878_30x_Q1 = 'data/NA12878/vcfeval/reaper.pacbio.blasr.30x.-z_-Q_1/{chrom}.done',
        NA12878_20x_Q50 = 'data/NA12878/vcfeval/reaper.pacbio.blasr.20x.-z_-Q_50/{chrom}.done',
        NA12878_20x_Q20 = 'data/NA12878/vcfeval/reaper.pacbio.blasr.20x.-z_-Q_20/{chrom}.done',
        NA12878_20x_Q1 = 'data/NA12878/vcfeval/reaper.pacbio.blasr.20x.-z_-Q_1/{chrom}.done',
    output:
        png = 'data/plots/effect_of_Q_param_PRcurve_{chrom}.png'
    run:
        ptf.plot_vcfeval([input.NA12878_30x_Q50[:-5], input.NA12878_20x_Q50[:-5],
                          input.NA12878_30x_Q20[:-5], input.NA12878_20x_Q20[:-5],
                          input.NA12878_30x_Q1[:-5], input.NA12878_20x_Q1[:-5]],
                         ['NA12878, 30x (Q50)', 'NA12878, 20x (Q50)',
                          'NA12878, 30x (Q20)', 'NA12878, 20x (Q20)',
                          'NA12878, 30x (Q1)', 'NA12878, 20x (Q1)'],
                           output.png,params.title,
                           colors=['#ff0000','#ff4f4f', # red
                                   '#ffff02','#ffff6b', #blue
                                   '#0400fc','#514efc'], #yellow
                           xlim=(0.85,1.0),ylim=(0.985,1.0))

rule plot_pr_curve_3_mappers:
    params: job_name = 'plot_pr_curve_3_mappers',
            title = 'Effect of Mapper on Precision/Recall for Reaper\n32x PacBio data for AJ Father, chr20'
    input:
        NA24149_blasr = 'data/NA24149/vcfeval/reaper.pacbio.blasr.32x.-z/20.done',
        NA24149_ngmlr = 'data/NA24149/vcfeval/reaper.pacbio.ngmlr.32x.-z/20.done',
        NA24149_bwa = 'data/NA24149/vcfeval/reaper.pacbio.bwa.32x.-z/20.done',
    output:
        png = 'data/plots/PR_curve_3_mappers_AJ_father_chr20.png'
    run:
        ptf.plot_vcfeval([input.NA24149_blasr[:-5],input.NA24149_ngmlr[:-5],input.NA24149_bwa[:-5]],
                         ['BLASR', 'NGMLR', 'BWA'],
                           output.png,params.title,
                           colors=['k','b','r'],
                           xlim=(0.5,1.0),ylim=(0.5,1.0))

rule plot_pr_curve_impact_of_haplotyping:
    params: job_name = 'plot_pr_curve_impact_of_haplotyping',
            title = 'Impact of Haplotype Information on PacBio Variant Calling with Reaper'
    input:
        NA12878_hap = 'data/NA12878/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA24385_hap = 'data/NA24385/vcfeval/reaper.pacbio.ngmlr.69x.-z/{chrom}.done',
        NA24149_hap = 'data/NA24149/vcfeval/reaper.pacbio.ngmlr.32x.-z/{chrom}.done',
        NA24143_hap = 'data/NA24143/vcfeval/reaper.pacbio.ngmlr.30x.-z/{chrom}.done',
        NA12878_no_hap = 'data/NA12878/vcfeval_no_haps/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA24385_no_hap = 'data/NA24385/vcfeval_no_haps/reaper.pacbio.ngmlr.69x.-z/{chrom}.done',
        NA24149_no_hap = 'data/NA24149/vcfeval_no_haps/reaper.pacbio.ngmlr.32x.-z/{chrom}.done',
        NA24143_no_hap = 'data/NA24143/vcfeval_no_haps/reaper.pacbio.ngmlr.30x.-z/{chrom}.done'
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

rule make_four_genomes_table:
    params: job_name = 'make_four_genomes_table.{chrom}.{GQ}'
    input:
        NA12878_vcfeval = 'data/NA12878/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA12878_vcfstats_genome = 'data/NA12878/variants/reaper.pacbio.blasr.44x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA12878_vcfstats_outside_giab = 'data/NA12878/variants/reaper.pacbio.blasr.44x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA12878_runtime = 'data/NA12878/variants/reaper.pacbio.blasr.44x.-z/{chrom}.vcf.runtime',
        NA24385_vcfeval = 'data/NA24385/vcfeval/reaper.pacbio.ngmlr.69x.-z/{chrom}.done',
        NA24385_vcfstats_genome = 'data/NA24385/variants/reaper.pacbio.ngmlr.69x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA24385_vcfstats_outside_giab = 'data/NA24385/variants/reaper.pacbio.ngmlr.69x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA24385_runtime = 'data/NA24385/variants/reaper.pacbio.ngmlr.69x.-z/{chrom}.vcf.runtime',
        NA24149_vcfeval = 'data/NA24149/vcfeval/reaper.pacbio.ngmlr.32x.-z/{chrom}.done',
        NA24149_vcfstats_genome = 'data/NA24149/variants/reaper.pacbio.ngmlr.32x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA24149_vcfstats_outside_giab = 'data/NA24149/variants/reaper.pacbio.ngmlr.32x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA24149_runtime = 'data/NA24149/variants/reaper.pacbio.ngmlr.32x.-z/{chrom}.vcf.runtime',
        NA24143_vcfeval = 'data/NA24143/vcfeval/reaper.pacbio.ngmlr.30x.-z/{chrom}.done',
        NA24143_vcfstats_genome = 'data/NA24143/variants/reaper.pacbio.ngmlr.30x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA24143_vcfstats_outside_giab = 'data/NA24143/variants/reaper.pacbio.ngmlr.30x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA24143_runtime = 'data/NA24143/variants/reaper.pacbio.ngmlr.30x.-z/{chrom}.vcf.runtime',
    output:
        table = 'data/output/four_GIAB_genomes_table.{chrom}.GQ{GQ}.tex'
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
                                 float(wildcards.GQ), output.table)

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
    output: vcfgz = 'data/{dataset}/variants/{calls_name}/{chrom}.GQ{GQ}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -g {wildcards.GQ} -i {input.vcfgz} -o {output.vcfgz}'

sim_covs = [20,30,40,60]
rule plot_simulation_pr_bars:
    params: job_name = 'plot_simulation_pr_bars.{chrom}.GQ{GQ}'
    input:
        reaper_genome = expand('data/simulation/vcfeval/reaper.pacbio.bwa.{c}x.-z/{{chrom}}.done',c=sim_covs),
        illumina_genome = expand('data/simulation/vcfeval/illumina_{c}x.filtered/{{chrom}}.done',c=sim_covs),
        reaper_segdup = expand('data/simulation/vcfeval_segdup/reaper.pacbio.bwa.{c}x.-z/{{chrom}}.done',c=sim_covs),
        illumina_segdup = expand('data/simulation/vcfeval_segdup/illumina_{c}x.filtered/{{chrom}}.done',c=sim_covs),
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
