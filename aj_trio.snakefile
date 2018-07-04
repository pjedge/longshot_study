rule aj_trio_mendelian_confident:
    params: job_name = 'aj_trio_mendelian.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz.tbi',
            ped = 'aj-trio.ped',
            sdf = 'data/genomes/1000g.sdf'
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/mendelian/{chrom}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g mendelian --lenient --pedigree={input.ped} -t {input.sdf} -i {input.vcfgz} --output={output.vcfgz}'

rule filter_trio_gq50_confident:
    params: job_name = 'filter_trio_gq50.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --all-samples -g 50 -i {input.vcfgz} -o {output.vcfgz}'

rule merge_SNVs_gt_20_trio_confident:
    params: job_name = 'merge_SNVs_confident_gt_20_trio.{chrom}'
    input:  vcfgz1 = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24385/reaper.pacbio.bwamem.69x.-z/{chrom}.vcf.gz',
            vcfgz2 = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24143/reaper.pacbio.bwamem.30x.-z/{chrom}.vcf.gz',
            vcfgz3 = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24149/reaper.pacbio.bwamem.32x.-z/{chrom}.vcf.gz',
            tbi1 = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24143/reaper.pacbio.bwamem.30x.-z/{chrom}.vcf.gz.tbi',
            tbi2 = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24149/reaper.pacbio.bwamem.32x.-z/{chrom}.vcf.gz.tbi',
            tbi3 = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24385/reaper.pacbio.bwamem.69x.-z/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz', #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.with_empties.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_gt_20_trio_confident:
    params: job_name = 'filter_SNVs_gt_20_trio.{id}.1000g.{callset}.{chrom}'
    input:  region_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/aj_trio_shared_confident_regions.all3_cov_greater_than_20.bed.gz',
            region_tbi = 'data/aj_trio/duplicated_regions/trio_covered_regions/aj_trio_shared_confident_regions.all3_cov_greater_than_20.bed.gz.tbi',
            vcfgz = 'data/{id}.1000g/variants/{callset}/{chrom}.vcf.gz',
            tbi = 'data/{id}.1000g/variants/{callset}/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/{id}/{callset}/{chrom, (\d+|all)}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.region_bed} -i {input.vcfgz} -o {output.vcfgz}'

rule get_aj_trio_cov_gt_20_confident:
    params: job_name = 'get_aj_trio_confident_cov_gt_20'
    input:  NA24143_bed = 'data/NA24143.1000g/variants/ground_truth/region_filter.bed',
            NA24149_bed = 'data/NA24149.1000g/variants/ground_truth/region_filter.bed',
            NA24385_bed = 'data/NA24385.1000g/variants/ground_truth/region_filter.bed',
            gt20_bed =  'data/aj_trio/duplicated_regions/trio_covered_regions/all.cov_greater_than_20.trio_intersect.bed'
    output: parents_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/aj_parents_shared_confident_regions.bed',
            trio_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/aj_trio_shared_confident_regions.bed',
            gt20_confident_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/aj_trio_shared_confident_regions.all3_cov_greater_than_20.bed',
    shell:
        '''
        {BEDTOOLS} intersect -sorted -a {input.NA24143_bed} -b {input.NA24149_bed} > {output.parents_bed}
        {BEDTOOLS} intersect -sorted -a {output.parents_bed} -b {input.NA24385_bed} > {output.trio_bed}
        {BEDTOOLS} intersect -sorted -a {output.trio_bed} -b {input.gt20_bed} > {output.gt20_confident_bed}
        '''

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule combine_aj_trio_cov_gt_20_all_confident:
    params: job_name = 'combine_aj_trio_cov_gt_20_all_confident'
    input:  expand('data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_20.trio_intersect.bed',chrom=chroms)
    output: 'data/aj_trio/duplicated_regions/trio_covered_regions/all.cov_greater_than_20.trio_intersect.bed'
    shell: 'rm -f {output}; for f in {input}; do {BEDTOOLS} sort -i $f >> {output}; done;'

#########################################################################################################################################
#########################################################################################################################################


# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule bgzip_trio_bed:
    params: job_name = 'bgzip_trio_bed'
    input:  'data/aj_trio/duplicated_regions/trio_covered_regions/{x}.bed'
    output: 'data/aj_trio/duplicated_regions/trio_covered_regions/{x}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule aj_trio_mendelian:
    params: job_name = 'aj_trio_mendelian.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz.tbi',
            ped = 'aj-trio.ped',
            sdf = 'data/genomes/1000g.sdf'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/{chrom}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g mendelian --lenient --pedigree={input.ped} -t {input.sdf} -i {input.vcfgz} --output={output.vcfgz}'

rule filter_trio_gq50:
    params: job_name = 'filter_trio_gq50.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --all-samples -g 50 -i {input.vcfgz} -o {output.vcfgz}'

rule merge_SNVs_gt_20_trio:
    params: job_name = 'merge_SNVs_gt_20_trio.{chrom}'
    input:  vcfgz1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.bwamem.69x.-z/{chrom}.vcf.gz',
            vcfgz2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.bwamem.30x.-z/{chrom}.vcf.gz',
            vcfgz3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.bwamem.32x.-z/{chrom}.vcf.gz',
            tbi1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.bwamem.30x.-z/{chrom}.vcf.gz.tbi',
            tbi2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.bwamem.32x.-z/{chrom}.vcf.gz.tbi',
            tbi3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.bwamem.69x.-z/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz', #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.with_empties.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_gt_20_trio:
    params: job_name = 'filter_SNVs_gt_20_trio.{id}.1000g.{callset}.{chrom}'
    input:  region_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_20.trio_intersect.segmental_duplications.bed.gz',
            vcfgz = 'data/{id}.1000g/variants/{callset}/{chrom}.vcf.gz',
            tbi = 'data/{id}.1000g/variants/{callset}/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id}/{callset}/{chrom, (\d+|all)}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.region_bed} -i {input.vcfgz} -o {output.vcfgz}'

rule combine_aj_trio_cov_gt_20:
    params: job_name = 'combine_aj_trio_cov_gt_20'
    input:  expand('data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_20.trio_intersect.segmental_duplications.bed',chrom=chroms)
    output: 'data/aj_trio/duplicated_regions/trio_covered_regions/all.cov_greater_than_20.trio_intersect.segmental_duplications.bed'
    shell: 'cat {input} > {output}'

rule get_aj_trio_cov_gt_20:
    params: job_name = 'get_aj_trio_cov_gt_20.{chrom}'
    input:  NA24143_bed = 'data/NA24143.1000g/aligned_reads/pacbio/pacbio.bwamem.{chrom}.30x.cov_greater_than_20.bed',
            NA24149_bed = 'data/NA24149.1000g/aligned_reads/pacbio/pacbio.bwamem.{chrom}.32x.cov_greater_than_20.bed',
            NA24385_bed = 'data/NA24385.1000g/aligned_reads/pacbio/pacbio.bwamem.{chrom}.69x.cov_greater_than_20.bed',
            segdup_bed =  'genome_tracks/segmental_duplications_0.99_similar_1000g.bed'
    output: parents_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom,\d+}.cov_greater_than_20.parents_intersect.bed',
            trio_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom,\d+}.cov_greater_than_20.trio_intersect.bed',
            trio_segdup_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom,\d+}.cov_greater_than_20.trio_intersect.segmental_duplications.bed',
    shell:
        '''
        {BEDTOOLS} intersect -a {input.NA24143_bed} -b {input.NA24149_bed} > {output.parents_bed}
        {BEDTOOLS} intersect -a {output.parents_bed} -b {input.NA24385_bed} > {output.trio_bed}
        {BEDTOOLS} intersect -a {output.trio_bed} -b {input.segdup_bed} > {output.trio_segdup_bed}
        '''

rule generate_coverage_bed:
    params: job_name = 'generate_coverage_bed.cov_greater_than_20.{id}.{build}.{chrom}.{cov}x'
    input:  'data/{id}.{build}/aligned_reads/pacbio/pacbio.bwamem.{chrom}.{cov}x.bam'
    output: 'data/{id}.{build}/aligned_reads/pacbio/pacbio.bwamem.{chrom,(\d+)}.{cov}x.cov_greater_than_20.bed',
    run:
        maxcov = int(wildcards.cov)*2 # exclude regions more than twice the mean coverage
        shell('''
        {SAMTOOLS} view -F 3844 -q 30 {input} -hb | \
        {BEDTOOLS} genomecov -bga -ibam - | \
        awk '$4 > 20 && $4 < {maxcov}' | \
        {BEDTOOLS} merge -i - > {output}
        ''')
