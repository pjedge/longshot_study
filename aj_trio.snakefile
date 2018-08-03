
from create_bed_poor_map import create_bed_poor_map

rule aj_trio_mendelian_confident:
    params: job_name = 'aj_trio_mendelian_confident.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz.tbi',
            ped = 'aj-trio.ped',
            sdf = 'data/genomes/hg38.sdf'
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/mendelian/{chrom}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g mendelian --lenient --pedigree={input.ped} -t {input.sdf} -i {input.vcfgz} --output={output.vcfgz}'

rule filter_trio_gq50_confident:
    params: job_name = 'filter_trio_gq50_confident.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.gq50.vcf.gz'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --all-samples -g 50 -i {input.vcfgz} -o {output.vcfgz}'

rule merge_SNVs_gt_20_trio_confident:
    params: job_name = 'merge_SNVs_confident_gt_20_trio'
    input:  vcfgz1 = expand('data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24385/reaper.pacbio.blasr.69x.-z/{chrom}.vcf.gz',chrom=chroms),
            vcfgz2 = expand('data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24143/reaper.pacbio.blasr.30x.-z/{chrom}.vcf.gz',chrom=chroms),
            vcfgz3 = expand('data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24149/reaper.pacbio.blasr.32x.-z/{chrom}.vcf.gz',chrom=chroms),
            tbi1 = expand('data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24143/reaper.pacbio.blasr.30x.-z/{chrom}.vcf.gz.tbi',chrom=chroms),
            tbi2 = expand('data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24149/reaper.pacbio.blasr.32x.-z/{chrom}.vcf.gz.tbi',chrom=chroms),
            tbi3 = expand('data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/NA24385/reaper.pacbio.blasr.69x.-z/{chrom}.vcf.gz.tbi',chrom=chroms)
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/merged/{chrom}.unfiltered.vcf.gz', #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.with_empties.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_gt_20_trio_confident:
    params: job_name = 'filter_SNVs_gt_20_trio.{id}.hg38.{callset}.{chrom}'
    input:  region_bed = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/aj_trio_shared_confident_regions.all3_cov_greater_than_20.chr{chrom}.bed.gz',
            region_tbi = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/aj_trio_shared_confident_regions.all3_cov_greater_than_20.chr{chrom}.bed.gz.tbi',
            vcfgz = 'data/{id}.hg38/variants/{callset}/{chrom}.vcf.gz',
            tbi = 'data/{id}.hg38/variants/{callset}/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/{id}/{callset}/{chrom, (\d+)}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=55g vcffilter --bed-regions={input.region_bed} -i {input.vcfgz} -o {output.vcfgz}'

rule get_aj_trio_cov_gt_20_confident:
    params: job_name = 'get_aj_trio_confident_cov_gt_20.{chrom}'
    input:  trio_bed = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/aj_trio_shared_confident_regions.bed',
            gt20_bed =  'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_20.trio_intersect.bed'
    output: gt20_confident_bed = 'data/aj_trio/duplicated_regions/test_confident_trio_shared_variant_sites/aj_trio_shared_confident_regions.all3_cov_greater_than_20.chr{chrom}.bed',
    shell:
        '''
        {BEDTOOLS} intersect -a {input.trio_bed} -b {input.gt20_bed} \
         | sort -k 1,1 -k2,2n > {output.gt20_confident_bed}
        '''

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
#rule combine_aj_trio_cov_gt_20_all_confident:
#    params: job_name = 'combine_aj_trio_cov_gt_20_all_confident'
#    input:  expand('data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_20.trio_intersect.bed',chrom=chroms)
#    output: 'data/aj_trio/duplicated_regions/trio_covered_regions/all.cov_greater_than_20.trio_intersect.bed'
#    shell: 'rm -f {output}; for f in {input}; do {BEDTOOLS} sort -i $f >> {output}; done;'

#########################################################################################################################################
#########################################################################################################################################


# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule bgzip_trio_bed:
    params: job_name = lambda wildcards: 'bgzip_trio_bed.{}'.format(str(wildcards.x).replace("/", "."))
    input:  'data/aj_trio/duplicated_regions/{x}.bed'
    output: 'data/aj_trio/duplicated_regions/{x}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule aj_trio_mendelian:
    params: job_name = 'aj_trio_mendelian.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.filtered.good_map_fraction.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.vcf.gz.tbi',
            ped = 'aj-trio.ped',
            sdf = 'data/genomes/hg38.sdf'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian/{chrom}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g mendelian --lenient --pedigree={input.ped} -t {input.sdf} -i {input.vcfgz} --output={output.vcfgz}'

rule filter_low_map_fraction_SNVs:
    params: job_name = 'filter_low_map_fraction_SNVs.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.filtered.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.vcf.gz.tbi',
            bed_file = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/merged.low_map_fraction_sites.bed'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.filtered.good_map_fraction.vcf.gz'
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-bed={input.bed_file} -i {input.vcfgz} -o {output.vcfgz}'

rule merge_low_map_frac_beds:
    params: job_name = 'merge_low_map_frac_beds'
    input: 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/father.low_map_fraction_sites.bed',
           'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/father.low_map_fraction_sites.bed',
    output: 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/merged.low_map_fraction_sites.bed',
    shell: 'cat {input} | {BEDTOOLS} sort -faidx genome_tracks/hg38.chrom.sizes.txt -i - | {BEDTOOLS} merge -i - > {output}'

rule determine_parents_low_map_fraction_sites:
    params: job_name = 'determine_parents_low_map_fraction_sites'
    input:  pileup_file_mother = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/NA24143.30x.pileup_over_mendelian_positions.txt',
            pileup_file_father = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/NA24149.32x.pileup_over_mendelian_positions.txt',
    output: bed_file_mother = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/mother.low_map_fraction_sites.bed',
            bed_file_father = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/mendelian.bak/father.low_map_fraction_sites.bed',
    run:
        create_bed_poor_map(pileup_file=input.pileup_file_mother,bed_file=output.bed_file_mother)
        create_bed_poor_map(pileup_file=input.pileup_file_father,bed_file=output.bed_file_father)

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_merged_SNVs:
    params: job_name = 'filter_merged_SNVs.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.vcf.gz.tbi',
            failed = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.merged_failed.vcf.gz'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.filtered.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-vcf={input.failed} -i {input.vcfgz} -o {output.vcfgz}'

rule merge_SNVs_gt_20_trio:
    params: job_name = 'merge_SNVs_gt_20_trio.{chrom}'
    input:  vcfgz1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.blasr.69x.-z/{chrom}.filtered.vcf.gz',
            vcfgz2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.blasr.30x.-z/{chrom}.filtered.vcf.gz',
            vcfgz3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.blasr.32x.-z/{chrom}.filtered.vcf.gz',
            tbi1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.blasr.30x.-z/{chrom}.filtered.vcf.gz.tbi',
            tbi2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.blasr.32x.-z/{chrom}.filtered.vcf.gz.tbi',
            tbi3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.blasr.69x.-z/{chrom}.filtered.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.vcf.gz', #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.with_empties.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'


rule merge_failed_SNVs:
    params: job_name = 'merge_failed_SNVs.{chrom}'
    input:  vcfgz1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.blasr.69x.-z/{chrom}.failed.vcf.gz',
            vcfgz2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.blasr.30x.-z/{chrom}.failed.vcf.gz',
            vcfgz3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.blasr.32x.-z/{chrom}.failed.vcf.gz',
            tbi1 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24143/reaper.pacbio.blasr.30x.-z/{chrom}.failed.vcf.gz.tbi',
            tbi2 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24149/reaper.pacbio.blasr.32x.-z/{chrom}.failed.vcf.gz.tbi',
            tbi3 = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/NA24385/reaper.pacbio.blasr.69x.-z/{chrom}.failed.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.merged_failed.vcf.gz', #'data/aj_trio/duplicated_regions/trio_shared_variant_sites/merged/{chrom}.unfiltered.with_empties.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'


# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_passing_SNVs:
    params: job_name = 'filter_passing_SNVs.{id}.hg38.{callset}.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id}/{callset}/{chrom}.filtered.vcf.gz',
            tbi   = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id}/{callset}/{chrom}.filtered.vcf.gz.tbi',
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id,(NA24143|NA24149|NA24385)}/{callset}/{chrom, (\d+|all)}.failed.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -r PASS -i {input.vcfgz} -o {output.vcfgz}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_trio:
    params: job_name = 'filter_SNVs_MEC_AQ_trio.{id}.hg38.{callset}.{chrom}'
    input:  vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id}/{callset}/{chrom}.vcf.gz',
            tbi = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id}/{callset}/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id,(NA24143|NA24149|NA24385)}/{callset}/{chrom, (\d+|all)}.filtered.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -g 50 --keep-expr "INFO.AQ > 7 && INFO.MB < 0.05 && INFO.MF < 0.1" --fail=fail -i {input.vcfgz} -o {output.vcfgz}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_gt_20_trio:
    params: job_name = 'filter_SNVs_gt_20_trio.{id}.hg38.{callset}.{chrom}'
    input:  region_bed = 'data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_20.trio_intersect.segmental_duplications.bed.gz',
            vcfgz = 'data/{id}.hg38/variants/{callset}/{chrom}.vcf.gz',
            tbi = 'data/{id}.hg38/variants/{callset}/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/duplicated_regions/trio_shared_variant_sites/{id}/{callset}/{chrom, (\d+|all)}.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.region_bed} -i {input.vcfgz} -o {output.vcfgz}'

rule combine_aj_trio_cov_gt_20:
    params: job_name = 'combine_aj_trio_cov_gt_20'
    input:  expand('data/aj_trio/duplicated_regions/trio_covered_regions/{chrom}.cov_greater_than_20.trio_intersect.segmental_duplications.bed',chrom=chroms)
    output: 'data/aj_trio/duplicated_regions/trio_covered_regions/all.cov_greater_than_20.trio_intersect.segmental_duplications.bed'
    shell: 'cat {input} > {output}'

rule get_aj_trio_cov_gt_20:
    params: job_name = 'get_aj_trio_cov_gt_20.{chrom}'
    input:  NA24143_bed = 'data/NA24143.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom}.30x.cov_greater_than_20.bed',
            NA24149_bed = 'data/NA24149.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom}.32x.cov_greater_than_20.bed',
            NA24385_bed = 'data/NA24385.hg38/aligned_reads/pacbio/pacbio.blasr.{chrom}.69x.cov_greater_than_20.bed',
            segdup_bed =  'genome_tracks/segmental_duplications_0.95_similar_hg38.bed'
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
    input:  'data/{id}.{build}/aligned_reads/pacbio/pacbio.blasr.{chrom}.{cov}x.bam'
    output: 'data/{id}.{build}/aligned_reads/pacbio/pacbio.blasr.{chrom,(\d+)}.{cov}x.cov_greater_than_20.bed',
    run:
        maxcov = int(wildcards.cov)*2 # exclude regions more than twice the mean coverage
        shell('''
        {SAMTOOLS} view -F 3844 -q 30 {input} -hb | \
        {BEDTOOLS} genomecov -bga -ibam - | \
        awk '$4 > 20 && $4 < {maxcov}' | \
        {BEDTOOLS} sort -i - | \
        {BEDTOOLS} merge -i - | \
        {BEDTOOLS} sort -i - > {output}
        ''')
