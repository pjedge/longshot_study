from paper_tables_and_figures import genomes_table_files

rule make_four_genomes_table:
    params: job_name = 'make_four_genomes_table.{chrom}.{GQ}'
    input:
        NA12878_vcfeval_dir = 'data/NA12878/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        NA12878_vcfstats_genome = 'data/NA12878/variants/reaper.pacbio.blasr.44x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA12878_vcfstats_outside_giab = 'data/NA12878/variants/reaper.pacbio.blasr.44x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA12878_runtime = 'data/NA12878/variants/reaper.pacbio.blasr.44x.-z/{chrom}.vcf.runtime',
        NA24385_vcfeval_dir = 'data/NA24385/vcfeval/reaper.pacbio.ngmlr.69x.-z/{chrom}.done',
        NA24385_vcfstats_genome = 'data/NA24385/variants/reaper.pacbio.ngmlr.69x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA24385_vcfstats_outside_giab = 'data/NA24385/variants/reaper.pacbio.ngmlr.69x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA24385_runtime = 'data/NA24385/variants/reaper.pacbio.ngmlr.69x.-z/{chrom}.vcf.runtime',
        NA24149_vcfeval_dir = 'data/NA24149/vcfeval/reaper.pacbio.ngmlr.32x.-z/{chrom}.done',
        NA24149_vcfstats_genome = 'data/NA24149/variants/reaper.pacbio.ngmlr.32x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA24149_vcfstats_outside_giab = 'data/NA24149/variants/reaper.pacbio.ngmlr.32x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA24149_runtime = 'data/NA24149/variants/reaper.pacbio.ngmlr.32x.-z/{chrom}.vcf.runtime',
        NA24143_vcfeval_dir = 'data/NA24143/vcfeval/reaper.pacbio.ngmlr.30x.-z/{chrom}.done',
        NA24143_vcfstats_genome = 'data/NA24143/variants/reaper.pacbio.ngmlr.30x.-z/{chrom}.GQ{GQ}.vcf.stats',
        NA24143_vcfstats_outside_giab = 'data/NA24143/variants/reaper.pacbio.ngmlr.30x.-z/{chrom}.outside_GIAB.GQ{GQ}.vcf.stats',
        NA24143_runtime = 'data/NA24143/variants/reaper.pacbio.ngmlr.30x.-z/{chrom}.vcf.runtime',
    output:
        table = 'data/output/four_GIAB_genomes_table.{chrom}.GQ{GQ}.tex'
    run:
        NA12878_table_files = genomes_table_files(input.NA12878_vcfeval_dir,
                                                  input.NA12878_vcfstats_genome,
                                                  input.NA12878_vcfstats_outside_giab,
                                                  input.NA12878_runtime)

        NA24385_table_files = genomes_table_files(input.NA24385_vcfeval_dir,
                                                  input.NA24385_vcfstats_genome,
                                                  input.NA24385_vcfstats_outside_giab,
                                                  input.NA24385_runtime)

        NA24149_table_files = genomes_table_files(input.NA24149_vcfeval_dir,
                                                  input.NA24149_vcfstats_genome,
                                                  input.NA24149_vcfstats_outside_giab,
                                                  input.NA24149_runtime)

        NA24143_table_files = genomes_table_files(input.NA24143_vcfeval_dir,
                                                  input.NA24143_vcfstats_genome,
                                                  input.NA24143_vcfstats_outside_giab,
                                                  input.NA24143_runtime)
        ptf.make_table_4_genomes(NA12878_table_files, NA24385_table_files,
                                 NA24149_table_files, NA24143_table_files,
                                 gq_cutoff=float(wildcards.GQ))

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
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-bed={input.region_filter} -g {wildcards.GQ} --no-index -i {input.vcfgz} -o {output.vcfgz}'

rule filter_vcf_GQ:
    params: job_name = 'filter_vcf_GQ.{dataset}.{calls_name}.{chrom}.GQ{GQ}'
    input:  vcfgz = 'data/{dataset}/variants/{calls_name}/{chrom}.vcf.gz'
    output: vcfgz = 'data/{dataset}/variants/{calls_name}/{chrom}.GQ{GQ}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -g {wildcards.GQ} --no-index -i {input.vcfgz} -o {output.vcfgz}'

sim_covs = [20,30,40,80]
rule plot_simulation_pr_bars:
    params: job_name = 'plot_simulation_pr_bars.{chrom}.GQ{GQ}'
    input:
        reaper_genome = expand('data/simulation/vcfeval/reaper.pacbio.bwa.{c}x.-z/{{chrom}}',c=sim_covs),
        illumina_genome = expand('data/simulation/vcfeval/illumina_{c}x.filtered/{{chrom}}',c=sim_covs),
        reaper_segdup = expand('data/simulation/vcfeval_segdup/reaper.pacbio.bwa.{c}x.-z/{{chrom}}',c=sim_covs),
        illumina_segdup = expand('data/simulation/vcfeval_segdup/illumina_{c}x.filtered/{{chrom}}',c=sim_covs),
    output:
        png = 'data/plots/simulation_pr_barplot_genome_vs_segdup.{chrom}.GQ{GQ}.png'
    run:
        ptf.plot_precision_recall_bars_simulation(
            list(input.reaper_genome),
            list(input.illumina_genome),
            list(input.reaper_segdup),
            list(input.illumina_segdup),
            float(wildcards.GQ),
            sim_covs,
            output.png)
