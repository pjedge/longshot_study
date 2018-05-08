#

rule plot_simulation_pr_bars:
    params: job_name = 'plot_simulation_pr_bars',
            title = 'Simulation PR'
    input:
        reaper30_rtg = 'data/NA12878/vcfeval/reaper.pacbio.blasr.30x.-z/{chrom}.done',
        reaper44_rtg = 'data/NA12878/vcfeval/reaper.pacbio.blasr.44x.-z/{chrom}.done',
        illumina_rtg = 'data/NA12878/vcfeval/illumina_30x.filtered/{chrom}.done'
    output:
        png = 'data/plots/NA12878_prec_recall_{chrom}.png'
    run:
        plot_vcfeval.plot_vcfeval(['data/NA12878/vcfeval/illumina_30x.filtered/{}'.format(wildcards.chrom),
                                   'data/NA12878/vcfeval/reaper.pacbio.blasr.30x.-z/{}'.format(wildcards.chrom),
                                   'data/NA12878/vcfeval/reaper.pacbio.blasr.44x.-z/{}'.format(wildcards.chrom)],
                                   ['Freebayes, Illumina 30x',
                                   'Reaper, PacBio 30x',
                                   'Reaper, PacBio 44x'],
                                   output.png,params.title,
                                   colors=['r','#8080ff','#3333ff'],
                                   xlim=(0.8,1.0),ylim=(0.985,1.0))
