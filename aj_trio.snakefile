
from create_bed_poor_map import create_bed_poor_map
import sys

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule bgzip_trio_bed:
    params: job_name = lambda wildcards: 'bgzip_trio_bed.{}.{}'.format(str(wildcards.aligner), str(wildcards.x).replace("/", "."))
    input:  'data/aj_trio.{aligner}/{x}.bed'
    output: 'data/aj_trio.{aligner}/{x}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule aj_trio_mendelian:
    params: job_name = 'aj_trio_mendelian.all.{tech}.{region}.{aligner}'
    input:  vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.filtered.no_centromeres.{region}.vcf.gz',
            tbi = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.filtered.no_centromeres.{region}.vcf.gz.tbi',
            ped = 'aj-trio.ped',
            sdf = 'data/genomes/hg38.sdf'
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/mendelian/all.{region}.vcf.gz',
            report = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/mendelian/all.{region}.report.txt',
    shell: '{RTGTOOLS} RTG_MEM=12g mendelian --lenient --pedigree={input.ped} -t {input.sdf} -i {input.vcfgz} --output={output.vcfgz} > {output.report}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_region_merged_SNVs:
    params: job_name = 'filter_region_merged_SNVs.all.{tech}.{region}.{aligner}'
    input:  vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.filtered.no_centromeres.vcf.gz',
            tbi = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.filtered.no_centromeres.vcf.gz.tbi',
            hg38_centromere_bed = 'genome_tracks/sv_repeat_telomere_centromere_hg38.bed',
            whole_genome_bed = 'genome_tracks/whole_genome_hg38.bed',
            confident_bed = 'data/aj_trio.{aligner}/confident_trio_intersect.bed',
            nonconfident_bed = 'data/aj_trio.{aligner}/nonconfident_trio_intersect.bed',
            segdup95_bed =  'genome_tracks/segmental_duplications_0.95_similar_hg38.bed',
            segdup99_bed =  'genome_tracks/segmental_duplications_0.99_similar_hg38.bed',
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/merged/all.filtered.no_centromeres.{region}.vcf.gz',
    run:
        if wildcards.region == 'whole_genome':
            region_bed = input.whole_genome_bed
        elif wildcards.region == 'confident':
            region_bed = input.confident_bed
        elif wildcards.region == 'nonconfident':
            region_bed = input.nonconfident_bed
        elif wildcards.region == 'segdup95':
            region_bed = input.segdup95_bed
        elif wildcards.region == 'segdup99':
            region_bed = input.segdup99_bed
        else:
            print("INVALID wildcard for \"region\"")
            exit(1)
        shell('{RTGTOOLS} RTG_MEM=12g vcffilter --include-bed={region_bed} -i {input.vcfgz} -o {output.vcfgz}')

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_centromere_merged_SNVs:
    params: job_name = 'filter_centromere_merged_SNVs.all.{tech}.{aligner}'
    input:  vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.filtered.vcf.gz',
            tbi = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.filtered.vcf.gz.tbi',
            hg38_centromere_bed = 'genome_tracks/sv_repeat_telomere_centromere_hg38.bed'
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/merged/all.filtered.no_centromeres.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-bed={input.hg38_centromere_bed} -i {input.vcfgz} -o {output.vcfgz}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_merged_SNVs:
    params: job_name = 'filter_failed_merged_SNVs.all.{tech}.{aligner}'
    input:  vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.vcf.gz',
            tbi = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.vcf.gz.tbi',
            failed = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/merged/all.merged_failed.vcf.gz',
            hg38_centromere_bed = 'genome_tracks/sv_repeat_telomere_centromere_hg38.bed'
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/merged/all.filtered.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-vcf={input.failed} -i {input.vcfgz} -o {output.vcfgz}'

rule merge_SNVs_gt_20_trio:
    params: job_name = 'merge_SNVs_gt_20_trio.all.{tech}.{aligner}'
    input:  vcfgz1 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24385/all.filtered.vcf.gz',
            vcfgz2 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24143/all.filtered.vcf.gz',
            vcfgz3 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24149/all.filtered.vcf.gz',
            tbi1 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24143/all.filtered.vcf.gz.tbi',
            tbi2 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24149/all.filtered.vcf.gz.tbi',
            tbi3 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24385/all.filtered.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/merged/all.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'


rule merge_failed_SNVs:
    params: job_name = 'merge_failed_SNVs.all.{tech}.{aligner}'
    input:  vcfgz1 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24385/all.failed.vcf.gz',
            vcfgz2 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24143/all.failed.vcf.gz',
            vcfgz3 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24149/all.failed.vcf.gz',
            tbi1 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24143/all.failed.vcf.gz.tbi',
            tbi2 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24149/all.failed.vcf.gz.tbi',
            tbi3 = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/NA24385/all.failed.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/merged/all.merged_failed.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_passing_SNVs:
    params: job_name = 'filter_passing_SNVs.hg38.{id}.{tech}.{aligner}.all'
    input:  vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/{id}/all.filtered.vcf.gz',
            tbi   = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech}/{id}/all.filtered.vcf.gz.tbi',
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(illumina|pacbio)}/{id,(NA24143|NA24149|NA24385)}/all.failed.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -r "PASS","." -i {input.vcfgz} -o {output.vcfgz}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_trio_illumina:
    params: job_name = 'filter_SNVs_MEC_AQ_trio_illumina.{id}.{aligner}.hg38.illumina.all'
    input:  vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/illumina/{id}/all.vcf.gz',
            tbi = 'data/aj_trio.{aligner}/trio_shared_variant_sites/illumina/{id}/all.vcf.gz.tbi'
    output: samplename_vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/illumina/{id,(NA24143|NA24149|NA24385)}/all.fixed_sample_name.vcf.gz',
            vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/illumina/{id,(NA24143|NA24149|NA24385)}/all.filtered.vcf.gz',
    shell:
        '''
        bcftools reheader -s <(echo "{wildcards.id}") {input.vcfgz} > {output.samplename_vcfgz}
        tabix -p vcf {output.samplename_vcfgz}
        {RTGTOOLS} RTG_MEM=12g vcffilter -g 50 --fail=fail -i {output.samplename_vcfgz} -o {output.vcfgz}
        '''

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_pacbio:
    params: job_name = 'filter_SNVs_MEC_AQ_trio_pacbio.{id}.{aligner}.hg38.pacbio.all'
    input:  vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/pacbio/{id}/all.vcf.gz',
            tbi = 'data/aj_trio.{aligner}/trio_shared_variant_sites/pacbio/{id}/all.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/pacbio/{id,(NA24143|NA24149|NA24385)}/all.filtered.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -g 50 --keep-expr "INFO.AQ > 7 && INFO.MB < 0.05 && INFO.MF < 0.1" --fail=fail -i {input.vcfgz} -o {output.vcfgz}'

callset_names = {
('pacbio','NA24385','minimap2') : 'longshot.pacbio.minimap2.69x._',
('pacbio','NA24149','minimap2') : 'longshot.pacbio.minimap2.32x._',
('pacbio','NA24143','minimap2') : 'longshot.pacbio.minimap2.30x._',
('illumina','NA24385','minimap2') : 'freebayes.illumina.aligned.60x.filtered',
('illumina','NA24149','minimap2') : 'freebayes.illumina.aligned.30x.filtered',
('illumina','NA24143','minimap2') : 'freebayes.illumina.aligned.30x.filtered',
('pacbio','NA24385','blasr') : 'longshot.pacbio.blasr.69x._',
('pacbio','NA24149','blasr') : 'longshot.pacbio.blasr.32x._',
('pacbio','NA24143','blasr') : 'longshot.pacbio.blasr.30x._',
('illumina','NA24385','blasr') : 'freebayes.illumina.aligned.60x.filtered',
('illumina','NA24149','blasr') : 'freebayes.illumina.aligned.30x.filtered',
('illumina','NA24143','blasr') : 'freebayes.illumina.aligned.30x.filtered'
}

dataset_aligners = {
('pacbio','blasr') : 'blasr',
('pacbio','minimap2') : 'minimap2',
('illumina','blasr') : 'aligned',       # we don't actually track aligner for short reads, they're just labeled as "aligned"
('illumina','minimap2') : 'aligned',
}

dataset_coverages = {
('pacbio','NA24385') : 69,
('pacbio','NA24149') : 32,
('pacbio','NA24143') : 30,
('illumina','NA24385') : 60,
('illumina','NA24149') : 30,
('illumina','NA24143') : 30
}

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_gt_20_trio:
    params: job_name = 'filter_SNVs_gt_20_trio.hg38.{id}.{tech}.{aligner}.all'
    input:  region_bed = 'data/aj_trio.{aligner}/trio_covered_regions/{tech}/all.cov_greater_than_20.trio_intersect.bed.gz',
            vcfgz = lambda wildcards: expand('data/{id}.hg38/variants/{x}/all.vcf.gz',id=wildcards.id,x=callset_names[(wildcards.tech,wildcards.id,wildcards.aligner)]),
            tbi = lambda wildcards: expand('data/{id}.hg38/variants/{x}/all.vcf.gz.tbi',id=wildcards.id,x=callset_names[(wildcards.tech,wildcards.id,wildcards.aligner)]),
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt'
    output: vcfgz = 'data/aj_trio.{aligner}/trio_shared_variant_sites/{tech,(pacbio|illumina)}/{id,(NA24385|NA24143|NA24149)}/all.vcf.gz',
    shell:
        '''
        {BEDTOOLS} intersect -u -wa -header -sorted -g {input.genome_file} -a {input.vcfgz} -b {input.region_bed} | \
        {BEDTOOLS} sort -g {input.genome_file} -header -i - | \
        bgzip -c > {output.vcfgz}
        '''
        #'{RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.region_bed} -i {input.vcfgz} -o {output.vcfgz}'

rule combine_aj_trio_cov_gt_20:
    params: job_name = 'combine_aj_trio_cov_gt_20.{aligner}.{tech}'
    input:  expand('data/aj_trio.{{aligner}}/trio_covered_regions/{{tech}}/{chrom}.cov_greater_than_20.trio_intersect.bed',chrom=chroms)
    output: 'data/aj_trio.{aligner}/trio_covered_regions/{tech,(pacbio|illumina)}/all.cov_greater_than_20.trio_intersect.bed'
    shell: 'cat {input} > {output}'

rule get_aj_trio_cov_gt_20:
    params: job_name = 'get_aj_trio_cov_gt_20.{aligner}.{tech}.{chrom}'
    input:  NA24143_bed = lambda wildcards: expand('data/NA24143.hg38/aligned_reads/{tech}/{tech}.{aligner}.{chrom}.{mother_cov}x.cov_greater_than_20.bed',tech=wildcards.tech, aligner=dataset_aligners[wildcards.tech, wildcards.aligner], chrom=wildcards.chrom, mother_cov=dataset_coverages[(wildcards.tech, 'NA24143')]),
            NA24149_bed = lambda wildcards: expand('data/NA24149.hg38/aligned_reads/{tech}/{tech}.{aligner}.{chrom}.{father_cov}x.cov_greater_than_20.bed',tech=wildcards.tech, aligner=dataset_aligners[wildcards.tech, wildcards.aligner], chrom=wildcards.chrom, father_cov=dataset_coverages[(wildcards.tech, 'NA24149')]),
            NA24385_bed = lambda wildcards: expand('data/NA24385.hg38/aligned_reads/{tech}/{tech}.{aligner}.{chrom}.{son_cov}x.cov_greater_than_20.bed',tech=wildcards.tech, aligner=dataset_aligners[wildcards.tech, wildcards.aligner], chrom=wildcards.chrom, son_cov=dataset_coverages[(wildcards.tech, 'NA24385')]),
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt' # not really needed here but eh
    output: parents_bed = 'data/aj_trio.{aligner}/trio_covered_regions/{tech,(pacbio|illumina)}/{chrom,\d+}.cov_greater_than_20.parents_intersect.bed',
            trio_bed = 'data/aj_trio.{aligner}/trio_covered_regions/{tech,(pacbio|illumina)}/{chrom,\d+}.cov_greater_than_20.trio_intersect.bed',
    shell:
        '''
        {BEDTOOLS} intersect -sorted -g {input.genome_file} -a {input.NA24143_bed} -b {input.NA24149_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.parents_bed}
        {BEDTOOLS} intersect -sorted -g {input.genome_file} -a {output.parents_bed} -b {input.NA24385_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.trio_bed}
        '''

rule merge_aj_trio_confident:
    params: job_name = 'merge_aj_trio_confident.{aligner}'
    input:  NA24143_bed = 'data/NA24143.hg38/variants/ground_truth/region_filter.bed',
            NA24149_bed = 'data/NA24149.hg38/variants/ground_truth/region_filter.bed',
            NA24385_bed = 'data/NA24385.hg38/variants/ground_truth/region_filter.bed',
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt'
    output: parents_bed = 'data/aj_trio.{aligner}/confident_parent_intersect.bed',
            trio_bed = 'data/aj_trio.{aligner}/confident_trio_intersect.bed',
    shell:
        '''
        {BEDTOOLS} intersect -a {input.NA24143_bed} -b {input.NA24149_bed} > {output.parents_bed}
        {BEDTOOLS} intersect -a {output.parents_bed} -b {input.NA24385_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.trio_bed}
        '''

rule merge_aj_trio_nonconfident:
    params: job_name = 'merge_aj_trio_nonconfident.{aligner}'
    input:  NA24143_bed = 'data/NA24143.hg38/variants/ground_truth/outside_region_filter.bed',
            NA24149_bed = 'data/NA24149.hg38/variants/ground_truth/outside_region_filter.bed',
            NA24385_bed = 'data/NA24385.hg38/variants/ground_truth/outside_region_filter.bed',
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt'
    output: parents_bed = 'data/aj_trio.{aligner}/nonconfident_parent_intersect.bed',
            trio_bed = 'data/aj_trio.{aligner}/nonconfident_trio_intersect.bed',
    shell:
        '''
        {BEDTOOLS} intersect -a {input.NA24143_bed} -b {input.NA24149_bed} > {output.parents_bed}
        {BEDTOOLS} intersect -a {output.parents_bed} -b {input.NA24385_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.trio_bed}
        '''

rule generate_coverage_bed:
    params: job_name = 'generate_coverage_bed_{tech}.{aligner}.cov_greater_than_20.{id}.{build}.{chrom}.{cov}x'
    input: bam = 'data/{id}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam',
           bai = 'data/{id}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam.bai',
           cov = 'data/{id}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{cov}x.bam.median_coverage'
    output: bed = 'data/{id}.{build}/aligned_reads/{tech}/{tech}.{aligner}.{chrom,(\d+)}.{cov}x.cov_greater_than_20.bed',
    run:
        median_cov = parse_int_file(input.cov)
        max_cov = int(median_cov + 4*sqrt(median_cov))

        shell('''
        {SAMTOOLS} view -F 3844 -q 30 {input.bam} chr{wildcards.chrom} -hb | \
        {BEDTOOLS} genomecov -bga -ibam - | \
        awk '$4 > 20 && $4 < {max_cov}' | \
        sort -k1,1 -k2,2n -S 1500M | \
        {BEDTOOLS} merge -i - | \
        sort -k1,1 -k2,2n -S 1500M > {output.bed}
        ''')
