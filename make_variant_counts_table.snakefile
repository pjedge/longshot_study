

rule make_variant_counts_table:
    params: job_name = 'make_variant_counts_table.{individual}.{build}.{aligner}.{pcov}.GQ{pGQ}',
    input: chr1_22_region_size_file = 'genome_tracks/whole_genome_{build}.bed.total_nonN__length',
           segmental_duplications_95_region_size_file = 'genome_tracks/segmental_duplications_0.95_similar_{build}.bed.total_nonN__length',
           segmental_duplications_99_region_size_file = 'genome_tracks/segmental_duplications_0.99_similar_{build}.bed.total_nonN__length',
           confident_region_size_file = 'data/{individual}.{build}/variants/ground_truth/region_filter.bed.total_nonN__length',
           nonconfident_region_size_file = 'data/{individual}.{build}/variants/ground_truth/outside_region_filter.bed.total_nonN__length',
           illumina_genome_stats = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.vcf.stats',
           illumina_segdup95_stats = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.segdup0.95_only.vcf.stats',
           illumina_segdup99_stats = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.segdup0.99_only.vcf.stats',
           illumina_GIAB_confident_stats = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.GIAB_confident_only.vcf.stats',
           illumina_GIAB_nonconfident_stats = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.stats',
           pacbio_genome_stats = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.stats',
           pacbio_segdup95_stats = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.segdup0.95_only.vcf.stats',
           pacbio_segdup99_stats = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.segdup0.99_only.vcf.stats',
           pacbio_GIAB_confident_stats = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.GIAB_confident_only.vcf.stats',
           pacbio_GIAB_nonconfident_stats = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.stats',
           pacbio_minus_illumina_genome_stats = 'data/{individual}.{build}/variants/MINUS_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}_freebayes.illumina.aligned.{icov}x.filtered_GQ50/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf.stats',
           pacbio_minus_illumina_segdup95_stats = 'data/{individual}.{build}/variants/MINUS_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}_freebayes.illumina.aligned.{icov}x.filtered_GQ50/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.segdup0.95_only.vcf.stats',
           pacbio_minus_illumina_segdup99_stats = 'data/{individual}.{build}/variants/MINUS_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}_freebayes.illumina.aligned.{icov}x.filtered_GQ50/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.segdup0.99_only.vcf.stats',
           pacbio_minus_illumina_GIAB_confident_stats = 'data/{individual}.{build}/variants/MINUS_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}_freebayes.illumina.aligned.{icov}x.filtered_GQ50/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_confident_only.vcf.stats',
           pacbio_minus_illumina_GIAB_nonconfident_stats = 'data/{individual}.{build}/variants/MINUS_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}_freebayes.illumina.aligned.{icov}x.filtered_GQ50/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.stats',
           illumina_minus_pacbio_genome_stats = 'data/{individual}.{build}/variants/MINUS_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf.stats',
           illumina_minus_pacbio_segdup95_stats = 'data/{individual}.{build}/variants/MINUS_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.segdup0.95_only.vcf.stats',
           illumina_minus_pacbio_segdup99_stats = 'data/{individual}.{build}/variants/MINUS_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.segdup0.99_only.vcf.stats',
           illumina_minus_pacbio_GIAB_confident_stats = 'data/{individual}.{build}/variants/MINUS_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_confident_only.vcf.stats',
           illumina_minus_pacbio_GIAB_nonconfident_stats = 'data/{individual}.{build}/variants/MINUS_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.stats',
           intersect_illumina_pacbio_genome_stats = 'data/{individual}.{build}/variants/INTERSECT_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/intersected.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf.stats',
           intersect_illumina_pacbio_segdup95_stats = 'data/{individual}.{build}/variants/INTERSECT_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/intersected.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.segdup0.95_only.vcf.stats',
           intersect_illumina_pacbio_segdup99_stats = 'data/{individual}.{build}/variants/INTERSECT_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/intersected.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.segdup0.99_only.vcf.stats',
           intersect_illumina_pacbio_GIAB_confident_stats = 'data/{individual}.{build}/variants/INTERSECT_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/intersected.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_confident_only.vcf.stats',
           intersect_illumina_pacbio_GIAB_nonconfident_stats = 'data/{individual}.{build}/variants/INTERSECT_freebayes.illumina.aligned.{icov}x.filtered_GQ50_longshot.pacbio.{aligner}.{pcov}x_GQ{pGQ}/intersected.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.stats',
    output: tex = 'data/output/variant_counts_table.{individual,(NA\d+)}.{build,(hg38|1000g)}.il{icov,\d+}x.{aligner,(blasr|minimap2)}.pb{pcov,\d+}x.GQ{pGQ,\d+}.tex',
    run:
        ptf.make_variant_counts_table(
               chr1_22_region_size_file = input.chr1_22_region_size_file,
               segmental_duplications_95_region_size_file = input.segmental_duplications_95_region_size_file,
               segmental_duplications_99_region_size_file = input.segmental_duplications_99_region_size_file,
               confident_region_size_file = input.confident_region_size_file,
               nonconfident_region_size_file = input.nonconfident_region_size_file,
               illumina_genome_stats = input.illumina_genome_stats,
               illumina_segdup95_stats = input.illumina_segdup95_stats,
               illumina_segdup99_stats = input.illumina_segdup99_stats,
               illumina_GIAB_confident_stats = input.illumina_GIAB_confident_stats,
               illumina_GIAB_nonconfident_stats = input.illumina_GIAB_nonconfident_stats,
               pacbio_genome_stats = input.pacbio_genome_stats,
               pacbio_segdup95_stats = input.pacbio_segdup95_stats,
               pacbio_segdup99_stats = input.pacbio_segdup99_stats,
               pacbio_GIAB_confident_stats = input.pacbio_GIAB_confident_stats,
               pacbio_GIAB_nonconfident_stats = input.pacbio_GIAB_nonconfident_stats,
               pacbio_minus_illumina_genome_stats = input.pacbio_minus_illumina_genome_stats,
               pacbio_minus_illumina_segdup95_stats = input.pacbio_minus_illumina_segdup95_stats,
               pacbio_minus_illumina_segdup99_stats = input.pacbio_minus_illumina_segdup99_stats,
               pacbio_minus_illumina_GIAB_confident_stats = input.pacbio_minus_illumina_GIAB_confident_stats,
               pacbio_minus_illumina_GIAB_nonconfident_stats = input.pacbio_minus_illumina_GIAB_nonconfident_stats,
               illumina_minus_pacbio_genome_stats = input.illumina_minus_pacbio_genome_stats,
               illumina_minus_pacbio_segdup95_stats = input.illumina_minus_pacbio_segdup95_stats,
               illumina_minus_pacbio_segdup99_stats = input.illumina_minus_pacbio_segdup99_stats,
               illumina_minus_pacbio_GIAB_confident_stats = input.illumina_minus_pacbio_GIAB_confident_stats,
               illumina_minus_pacbio_GIAB_nonconfident_stats = input.illumina_minus_pacbio_GIAB_nonconfident_stats,
               intersect_illumina_pacbio_genome_stats = input.intersect_illumina_pacbio_genome_stats,
               intersect_illumina_pacbio_segdup95_stats = input.intersect_illumina_pacbio_segdup95_stats,
               intersect_illumina_pacbio_segdup99_stats = input.intersect_illumina_pacbio_segdup99_stats,
               intersect_illumina_pacbio_GIAB_confident_stats = input.intersect_illumina_pacbio_GIAB_confident_stats,
               intersect_illumina_pacbio_GIAB_nonconfident_stats = input.intersect_illumina_pacbio_GIAB_nonconfident_stats,
               outfile=output.tex)

rule generate_nonconfident_bed:
    params: job_name = 'generate_nonconfident_bed.{individual}.{build}',
    input:  whole_genome_bed = 'genome_tracks/whole_genome_{build}.bed',
            confident_bed = 'data/{individual}.{build}/variants/ground_truth/region_filter.bed',
    output: bed = 'data/{individual}.{build}/variants/ground_truth/outside_region_filter.bed',
    shell: '{BEDTOOLS} subtract -a {input.whole_genome_bed} -b {input.confident_bed} > {output.bed}'

##################################################################################################
# IMPORTANT NOTE
# the next 6 rules use the -sorted option for faster processing in bedtools
# since we're operating on VCF files in natural 1..22 chromosome order,
# we need to make sure to pass a genome file that is also in this order
##################################################################################################

rule intersect_illumina_pacbio_VCFs:
    params: job_name = 'intersect_VCFs.{individual}.{build}.{aligner}.{pcov}x.illumina_GQ{iGQ}.pacbio_GQ{pGQ}',
    input:  pacbio_vcfgz =   'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            pacbio_tbi   =   'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            illumina_vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ{iGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            illumina_tbi =   'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ{iGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            genome_file = 'genome_tracks/{build}.chrom.sizes.natural_order.txt'
    output: vcf = 'data/{individual}.{build}/variants/INTERSECT_freebayes.illumina.aligned.{icov,\d+}x.filtered_GQ{iGQ,\d+}_longshot.pacbio.{aligner,(blasr|minimap2)}.{pcov,\d+}x_GQ{pGQ}/intersected.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf',
            vcfgz = 'data/{individual}.{build}/variants/INTERSECT_freebayes.illumina.aligned.{icov,\d+}x.filtered_GQ{iGQ,\d+}_longshot.pacbio.{aligner,(blasr|minimap2)}.{pcov,\d+}x_GQ{pGQ}/intersected.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
    shell:
        '''
        gunzip -c {input.pacbio_vcfgz} | grep -P '^#' > {output.vcf}
        {BEDTOOLS} intersect -g {input.genome_file} -sorted -wa -a {input.pacbio_vcfgz} -b {input.illumina_vcfgz} >> {output.vcf}
        bgzip -c {output.vcf} > {output.vcfgz}
        '''

rule setminus_pacbio_illumina_VCFs:
    params: job_name = 'setminus_pacbio_illumina_VCFs.{individual}.{build}.{aligner}.{pcov}x.illumina_GQ50.pacbio_GQ{pGQ}',
    input:  pacbio_vcfgz =   'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            pacbio_tbi   =   'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            illumina_vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            illumina_tbi =   'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            genome_file = 'genome_tracks/{build}.chrom.sizes.natural_order.txt'
    output: vcf = 'data/{individual}.{build}/variants/MINUS_longshot.pacbio.{aligner,(blasr|minimap2)}.{pcov,\d+}x_GQ{pGQ,\d+}_freebayes.illumina.aligned.{icov,\d+}x.filtered_GQ50/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf',
            vcfgz = 'data/{individual}.{build}/variants/MINUS_longshot.pacbio.{aligner,(blasr|minimap2)}.{pcov,\d+}x_GQ{pGQ,\d+}_freebayes.illumina.aligned.{icov,\d+}x.filtered_GQ50/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
    shell:
        '''
        gunzip -c {input.pacbio_vcfgz} | grep -P '^#' > {output.vcf}
        {BEDTOOLS} subtract -g {input.genome_file} -sorted -a {input.pacbio_vcfgz} -b {input.illumina_vcfgz} >> {output.vcf}
        bgzip -c {output.vcf} > {output.vcfgz}
        '''

rule setminus_illumina_pacbio_VCFs:
    params: job_name = 'setminus_illumina_pacbio_VCFs.{individual}.{build}.{aligner}.{pcov}x.illumina_GQ50.pacbio_GQ{pGQ}',
    input:  pacbio_vcfgz =   'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            pacbio_tbi   =   'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/all.GQ{pGQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            illumina_vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            illumina_tbi =   'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/all.GQ50.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            genome_file = 'genome_tracks/{build}.chrom.sizes.natural_order.txt'
    output: vcf = 'data/{individual}.{build}/variants/MINUS_freebayes.illumina.aligned.{icov,\d+}x.filtered_GQ50_longshot.pacbio.{aligner,(blasr|minimap2)}.{pcov,\d+}x_GQ{pGQ,\d+}/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf',
            vcfgz = 'data/{individual}.{build}/variants/MINUS_freebayes.illumina.aligned.{icov,\d+}x.filtered_GQ50_longshot.pacbio.{aligner,(blasr|minimap2)}.{pcov,\d+}x_GQ{pGQ,\d+}/minus.all.GQ0.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
    shell:
        '''
        gunzip -c {input.illumina_vcfgz} | grep -P '^#' > {output.vcf}
        {BEDTOOLS} subtract -g {input.genome_file} -sorted -a {input.illumina_vcfgz} -b {input.pacbio_vcfgz} >> {output.vcf}
        bgzip -c {output.vcf} > {output.vcfgz}
        '''

rule count_bed_size_genometrack:
    params: job_name = 'count_bed_size_genometrack.{x}_{build}'
    input: bed = 'genome_tracks/{x}_{build}.bed',
           nonN = 'genome_tracks/nonN_regions_{build}.bed'
    output: len = 'genome_tracks/{x}_{build,(hg38|1000g)}.bed.total_nonN__length'
    shell:
        '''
        {BEDTOOLS} intersect -a {input.bed} -b {input.nonN} | \
        {BEDTOOLS} sort -i - | \
        {BEDTOOLS} merge -i - | \
        python filter_bed_chroms.py | \
        awk -F'\t' 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}' \
        > {output}
        '''

rule count_bed_size_conf_nonconf:
    params: job_name = 'count_bed_size_conf_nonconf.{individual}.{build}.{x}'
    input: bed = 'data/{individual}.{build}/variants/ground_truth/{x}.bed',
           nonN = 'genome_tracks/nonN_regions_{build}.bed'
    output: len = 'data/{individual}.{build}/variants/ground_truth/{x}.bed.total_nonN__length'
    shell:
        '''
        {BEDTOOLS} intersect -a {input.bed} -b {input.nonN} | \
        {BEDTOOLS} sort -i - | \
        {BEDTOOLS} merge -i - | \
        python filter_bed_chroms.py | \
        awk -F'\t' 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}' \
        > {output}
        '''

rule get_nonN_bedfile:
    params: job_name = 'get_nonN_bedfile.{build}'
    input: n_regions = 'genome_tracks/N_regions_{build}.bed',
           whole_genome = 'genome_tracks/whole_genome_{build}.bed'
    output: nonN_regions = 'genome_tracks/nonN_regions_{build}.bed'
    shell: '{BEDTOOLS} subtract -a {input.whole_genome} -b {input.n_regions} > {output.nonN_regions}'

rule get_N_bedfile:
    params: job_name = 'get_N_bedfile.{build}'
    input: 'data/genomes/{build}.fa'
    output: 'genome_tracks/N_regions_{build}.bed'
    shell: '{SEQTK} cutN -gp10000000 -n1 {input} > {output}'

rule filter_SNVs_GQ:
    params: job_name = 'filter_SNVs_GQ.{individual}.{build}.{info}.GQ{GQ}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{info}/all.DECOMPOSED.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/{info}/all.DECOMPOSED.vcf.gz.tbi',
    output: vcfgz = 'data/{individual}.{build}/variants/{info}/all.GQ{GQ,\d+}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter \
        --keep-expr "FILTER.every(function(f) {{return f == 'PASS'}})" \
        --snps-only -g {wildcards.GQ} -i {input.vcfgz} -o {output.vcfgz}
        '''

rule decompose_SNVs:
    params: job_name = 'decompose_SNVs.{individual}.{build}.{info}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{info}/all.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/{info}/all.vcf.gz.tbi',
    output: vcfgz =  'data/{individual}.{build}/variants/{info}/all.DECOMPOSED.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfdecompose --break-mnps --break-indels -i {input.vcfgz} -o {output.vcfgz}'

rule filter_SNVs_segmental_duplications:
    params: job_name = 'filter_SNVs_segmental_duplications.{individual}.{build}.{info}.GQ{GQ}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            bed   = 'genome_tracks/segmental_duplications_{frac}_similar_{build}.bed.gz'
    output: vcfgz = 'data/{individual}.{build}/variants/{info}/{type,(all|intersected\.all|minus\.all)}.GQ{GQ,\d+}.PASS.SNPs_ONLY.DECOMPOSED.segdup{frac}_only.vcf.gz'
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.bed} \
        -i {input.vcfgz} -o {output.vcfgz}
        '''

rule filter_SNVs_coding_exons_segmental_duplications:
    params: job_name = 'filter_SNVs_coding_exons_segmental_duplications.{individual}.{build}.{info}.GQ{GQ}.{frac}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{info}/all.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/{info}/all.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            bed   = 'genome_tracks/coding_exons_intersect_segmental_duplications_{build}_{frac}_similar.bed.gz'
    output: vcfgz = 'data/{individual}.{build}/variants/{info}/all.GQ{GQ,\d+}.PASS.SNPs_ONLY.DECOMPOSED.coding_exons_segdup{frac}_only.vcf.gz'
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.bed} \
        -i {input.vcfgz} -o {output.vcfgz}
        '''

rule filter_SNVs_confident_regions:
    params: job_name = 'filter_SNVs_GIAB_confident.{individual}.{build}.{info}.GQ{GQ}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            bed   = 'data/{individual}.{build}/variants/ground_truth/region_filter.bed.gz'
    output: vcfgz = 'data/{individual}.{build}/variants/{info}/{type,(all|intersected\.all|minus\.all)}.GQ{GQ,\d+}.PASS.SNPs_ONLY.DECOMPOSED.GIAB_confident_only.vcf.gz'
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.bed} \
        -i {input.vcfgz} -o {output.vcfgz}
        '''

rule filter_SNVs_outside_confident_regions:
    params: job_name = 'filter_SNVs_outside_GIAB_confident.{individual}.{build}.{info}.GQ{GQ}',
    input:  vcfgz = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/{info}/{type}.GQ{GQ}.PASS.SNPs_ONLY.DECOMPOSED.vcf.gz.tbi',
            bed   = 'data/{individual}.{build}/variants/ground_truth/region_filter.bed.gz'
    output: vcfgz = 'data/{individual}.{build}/variants/{info}/{type,(all|intersected\.all|minus\.all)}.GQ{GQ,\d+}.PASS.SNPs_ONLY.DECOMPOSED.GIAB_nonconfident_only.vcf.gz'
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter --exclude-bed={input.bed} \
        -i {input.vcfgz} -o {output.vcfgz}
        '''
