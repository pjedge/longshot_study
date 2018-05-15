# filter variants for aj trio individuals for positions that have at least 30x coverage in all 3 datasets
rule aj_trio_mendelian:
    params: job_name = 'aj_trio_mendelian.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz.tbi',
            ped = 'aj-trio.ped',
            sdf = 'data/genomes/1000g_v37_phase2.sdf'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/{chrom}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g mendelian --lenient --pedigree={input.ped} -t {input.sdf} -i {input.vcfgz} --output={output.vcfgz}'

# filter variants for aj trio individuals for positions that have at least 30x coverage in all 3 datasets
rule filter_trio_gq50:
    params: job_name = 'filter_trio_gq50.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --all-samples -g 50 -i {input.vcfgz} -o {output.vcfgz}'

#rule replace_empty_gt_with_reference:
#    params: job_name = 'replace_empty_gt_with_reference.{chrom}'
#    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.with_empties.vcf.gz'
#    output: vcf = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom, (\d+|X|Y)}.unfiltered.vcf',
#            vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom, (\d+|X|Y)}.unfiltered.vcf.gz'
#    run:
#        replace_empty_gt_with_reference(input.vcfgz, output.vcf)
#        shell('{BGZIP} -c {output.vcf} > {output.vcfgz}')

# filter variants for aj trio individuals for positions that have at least 30x coverage in all 3 datasets
rule merge_SNVs_gt_30_trio:
    params: job_name = 'merge_SNVs_gt_30_trio.{chrom}'
    input:  vcfgz1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.ngmlr.69x.-z/{chrom}.vcf.gz',
            vcfgz2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.ngmlr.30x.-z/{chrom}.vcf.gz',
            vcfgz3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.ngmlr.32x.-z/{chrom}.vcf.gz',
            tbi1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.ngmlr.30x.-z/{chrom}.vcf.gz.tbi',
            tbi2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.ngmlr.32x.-z/{chrom}.vcf.gz.tbi',
            tbi3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.ngmlr.69x.-z/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz', #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.with_empties.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'

# filter variants for aj trio individuals for positions that have at least 30x coverage in all 3 datasets
rule filter_SNVs_gt_30_trio:
    params: job_name = 'filter_SNVs_gt_30_trio.{id}.{callset}.{chrom}'
    input:  region_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_30.trio_intersect.segmental_duplications.bed',
            vcfgz = 'data/{id}/variants/{callset}/{chrom}.vcf.gz',
            tbi = 'data/{id}/variants/{callset}/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id}/{callset}/{chrom, (\d+|X|Y)}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.region_bed} -i {input.vcfgz} -o {output.vcfgz}'

rule get_aj_trio_cov_gt_30:
    params: job_name = 'get_aj_trio_cov_gt_30.{chrom}'
    input:  NA24143_bed = 'data/NA24143/aligned_reads/pacbio/pacbio.ngmlr.{chrom}.30x.cov_greater_than_30.bed',
            NA24149_bed = 'data/NA24149/aligned_reads/pacbio/pacbio.ngmlr.{chrom}.32x.cov_greater_than_30.bed',
            NA24385_bed = 'data/NA24385/aligned_reads/pacbio/pacbio.ngmlr.{chrom}.69x.cov_greater_than_30.bed',
            segdup_bed =  'genome_tracks/segmental_duplications_0.99_similar_1000g.bed'
    output: parents_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_30.parents_intersect.bed',
            trio_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_30.trio_intersect.bed',
            trio_segdup_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_30.trio_intersect.segmental_duplications.bed',
    shell:
        '''
        {BEDTOOLS} intersect -a {input.NA24143_bed} -b {input.NA24149_bed} > {output.parents_bed}
        {BEDTOOLS} intersect -a {output.parents_bed} -b {input.NA24385_bed} > {output.trio_bed}
        {BEDTOOLS} intersect -a {output.trio_bed} -b {input.segdup_bed} > {output.trio_segdup_bed}
        '''

rule generate_coverage_bed:
    params: job_name = 'generate_coverage_bed.cov_greater_than_30.{id}.{chrom}.{cov}x'
    input:  'data/{id}/aligned_reads/pacbio/pacbio.ngmlr.{chrom}.{cov}x.bam'
    output: 'data/{id}/aligned_reads/pacbio/pacbio.ngmlr.{chrom}.{cov}x.cov_greater_than_30.bed',
    run:
        maxcov = int(wildcards.cov)*2
        shell('''
        {SAMTOOLS} view -F 3844 -q 30 {input} -hb | \
        {BEDTOOLS} genomecov -bga -ibam - | \
        awk '$4 > 30 && $4 < {maxcov}' | \
        {BEDTOOLS} merge -i - > {output}
        ''')
