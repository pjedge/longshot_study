
from create_bed_poor_map import create_bed_poor_map
import sys

rule generate_aj_trio_table:
    params: job_name = 'generate_aj_trio_table'
    input: expand('data/aj_trio/{region}.m30.f30.s60/trio_shared_variant_sites/mendelian/illumina/all.report.txt',region=['whole_genome','confident','nonconfident','segdup95','segdup99']),
           expand('data/aj_trio/{region}.m30.f32.s69/trio_shared_variant_sites/mendelian/pacbio/all.report.txt',region=['whole_genome','confident','nonconfident','segdup95','segdup99'])


# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule bgzip_trio_bed:
    params: job_name = lambda wildcards: 'bgzip_trio_bed.{}.m{}.f{}.s{}.{}'.format(wildcards.region, wildcards.m, wildcards.f, wildcards.s, str(wildcards.x).replace("/", "."))
    input:  'data/aj_trio/{region}.m{m}.f{f}.s{s}/{x}.bed'
    output: 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/{x}.bed.gz'
    shell: 'bgzip -c {input} > {output}'

rule aj_trio_mendelian:
    params: job_name = 'aj_trio_mendelian.all.{region}.m{m}.f{f}.s{s}'
    input:  vcfgz = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/merged/{tech}/all.filtered.vcf.gz',
            tbi = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/merged/{tech}/all.filtered.vcf.gz.tbi',
            ped = 'aj-trio.ped',
            sdf = 'data/genomes/hg38.sdf'
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/mendelian/{tech,(illumina|pacbio)}/all.vcf.gz',
            report = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/mendelian/{tech,(illumina|pacbio)}/all.report.txt',
    shell: '{RTGTOOLS} RTG_MEM=12g mendelian --lenient --pedigree={input.ped} -t {input.sdf} -i {input.vcfgz} --output={output.vcfgz} > {output.report}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_merged_SNVs:
    params: job_name = 'filter_merged_SNVs.all.{region}.m{m}.f{f}.s{s}'
    input:  vcfgz = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/merged/{tech}/all.vcf.gz',
            tbi = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/merged/{tech}/all.vcf.gz.tbi',
            failed = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/merged/{tech}/all.merged_failed.vcf.gz'
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/merged/{tech,(illumina|pacbio)}/all.filtered.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter --exclude-vcf={input.failed} -i {input.vcfgz} -o {output.vcfgz}'

rule merge_SNVs_gt_20_trio:
    params: job_name = 'merge_SNVs_gt_20_trio.all.{region}.m{m}.f{f}.s{s}'
    input:  vcfgz1 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24385/{tech}/all.filtered.vcf.gz',
            vcfgz2 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24143/{tech}/all.filtered.vcf.gz',
            vcfgz3 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24149/{tech}/all.filtered.vcf.gz',
            tbi1 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24143/{tech}/all.filtered.vcf.gz.tbi',
            tbi2 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24149/{tech}/all.filtered.vcf.gz.tbi',
            tbi3 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24385/{tech}/all.filtered.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/merged/{tech,(illumina|pacbio)}/all.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'


rule merge_failed_SNVs:
    params: job_name = 'merge_failed_SNVs.all.{region}.m{m}.f{f}.s{s}'
    input:  vcfgz1 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24385/{tech}/all.failed.vcf.gz',
            vcfgz2 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24143/{tech}/all.failed.vcf.gz',
            vcfgz3 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24149/{tech}/all.failed.vcf.gz',
            tbi1 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24143/{tech}/all.failed.vcf.gz.tbi',
            tbi2 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24149/{tech}/all.failed.vcf.gz.tbi',
            tbi3 = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/NA24385/{tech}/all.failed.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/merged/{tech,(illumina|pacbio)}/all.merged_failed.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcfmerge -o {output} {input.vcfgz1} {input.vcfgz2} {input.vcfgz3}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_passing_SNVs:
    params: job_name = 'filter_passing_SNVs.{id}.hg38.{region}.m{m}.f{f}.s{s}.{tech}.all'
    input:  vcfgz = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/{id}/{tech}/all.filtered.vcf.gz',
            tbi   = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/{id}/{tech}/all.filtered.vcf.gz.tbi',
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/{id,(NA24143|NA24149|NA24385)}/{tech,(illumina|pacbio)}/all.failed.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -r "PASS","." -i {input.vcfgz} -o {output.vcfgz}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_trio_illumina:
    params: job_name = 'filter_SNVs_MEC_AQ_trio_illumina.{id}.hg38.{region}.m{m}.f{f}.s{s}.illumina.all'
    input:  vcfgz = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/{id}/illumina/all.vcf.gz',
            tbi = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/{id}/illumina/all.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/{id,(NA24143|NA24149|NA24385)}/illumina/all.filtered.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -g 50 --fail=fail -i {input.vcfgz} -o {output.vcfgz}'

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_pacbio:
    params: job_name = 'filter_SNVs_MEC_AQ_trio_pacbio.{id}.hg38.{region}.m{m}.f{f}.s{s}.pacbio.all'
    input:  vcfgz = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/{id}/pacbio/all.vcf.gz',
            tbi = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_shared_variant_sites/{id}/pacbio/all.vcf.gz.tbi'
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/{id,(NA24143|NA24149|NA24385)}/pacbio/all.filtered.vcf.gz',
    shell: '{RTGTOOLS} RTG_MEM=12g vcffilter -g 50 --keep-expr "INFO.AQ > 7 && INFO.MB < 0.05 && INFO.MF < 0.1" --fail=fail -i {input.vcfgz} -o {output.vcfgz}'


def get_callset_name(id, tech, mother_cov, father_cov, son_cov):
    if tech == 'pacbio':
        if id == 'NA24385':
            return 'reaper.pacbio.blasr.{}x.-z'.format(son_cov)
        elif id == 'NA24149':
            return 'reaper.pacbio.blasr.{}x.-z'.format(father_cov)
        elif id == 'NA24143':
            return 'reaper.pacbio.blasr.{}x.-z'.format(mother_cov)
    elif tech == 'illumina':
        if id == 'NA24385':
            return 'illumina_{}x.filtered'.format(son_cov)
        elif id == 'NA24149':
            return 'illumina_{}x.filtered'.format(father_cov)
        elif id == 'NA24143':
            return 'illumina_{}x.filtered'.format(mother_cov)
    print("ERROR: unable to determine callset name for combination of individual, coverage, and technology:", file=sys.stderr)
    print("individual: {}, technology: {}".format(id, tech), file=sys.stderr)
    exit(1)

# filter variants for aj trio individuals for positions that have at least 20x coverage in all 3 datasets
rule filter_SNVs_gt_20_trio:
    params: job_name = 'filter_SNVs_gt_20_trio.{id}.hg38.{tech}.all.{region}.m{m}.f{f}.s{s}'
    input:  region_bed = 'data/aj_trio/{region}.m{m}.f{f}.s{s}/trio_covered_regions/{tech}/all.cov_greater_than_20.trio_intersect.region_filtered.bed.gz',
            vcfgz = lambda wildcards: expand('data/{id}.hg38/variants/{x}/all.vcf.gz',id=wildcards.id,x=get_callset_name(wildcards.id, wildcards.tech, wildcards.m, wildcards.f, wildcards.s)),
            tbi = lambda wildcards: expand('data/{id}.hg38/variants/{x}/all.vcf.gz.tbi',id=wildcards.id,x=get_callset_name(wildcards.id, wildcards.tech, wildcards.m, wildcards.f, wildcards.s)),
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt'
    output: vcfgz = 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_shared_variant_sites/{id,(NA24385|NA24143|NA24149)}/{tech,(pacbio|illumina)}/all.vcf.gz',
    shell:
        '''
        {BEDTOOLS} intersect -u -wa -header -sorted -g {input.genome_file} -a {input.vcfgz} -b {input.region_bed} | \
        {BEDTOOLS} sort -g {input.genome_file} -header -i - | \
        bgzip -c > {output.vcfgz}
        '''
        #'{RTGTOOLS} RTG_MEM=12g vcffilter --bed-regions={input.region_bed} -i {input.vcfgz} -o {output.vcfgz}'

rule combine_aj_trio_cov_gt_20:
    params: job_name = 'combine_aj_trio_cov_gt_20.{tech}.{region}.m{m}.f{f}.s{s}'
    input:  expand('data/aj_trio/{{region}}.m{{m}}.f{{f}}.s{{s}}/trio_covered_regions/{{tech}}/{chrom}.cov_greater_than_20.trio_intersect.region_filtered.bed',chrom=chroms)
    output: 'data/aj_trio/{region}.m{m,\d+}.f{f,\d+}.s{s,\d+}/trio_covered_regions/{tech,(pacbio|illumina)}/all.cov_greater_than_20.trio_intersect.region_filtered.bed'
    shell: 'cat {input} > {output}'

rule get_aj_trio_cov_gt_20:
    params: job_name = 'get_aj_trio_cov_gt_20.{tech}.{region}.{chrom}'
    input:  NA24143_bed = 'data/NA24143.hg38/aligned_reads/{tech}/{tech}.{chrom}.{mother_cov}x.cov_greater_than_20.bed',
            NA24149_bed = 'data/NA24149.hg38/aligned_reads/{tech}/{tech}.{chrom}.{father_cov}x.cov_greater_than_20.bed',
            NA24385_bed = 'data/NA24385.hg38/aligned_reads/{tech}/{tech}.{chrom}.{son_cov}x.cov_greater_than_20.bed',
            whole_genome_bed = 'genome_tracks/whole_genome_hg38.bed',
            confident_bed = 'data/aj_trio/confident_trio_intersect.bed',
            nonconfident_bed = 'data/aj_trio/nonconfident_trio_intersect.bed',
            segdup95_bed =  'genome_tracks/segmental_duplications_0.95_similar_hg38.bed',
            segdup99_bed =  'genome_tracks/segmental_duplications_0.99_similar_hg38.bed',
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt' # not really needed here but eh
    output: parents_bed = 'data/aj_trio/{region}.m{mother_cov}.f{father_cov}.s{son_cov}/trio_covered_regions/{tech,(pacbio|illumina)}/{chrom,\d+}.cov_greater_than_20.parents_intersect.bed',
            trio_bed = 'data/aj_trio/{region}.m{mother_cov}.f{father_cov}.s{son_cov}/trio_covered_regions/{tech,(pacbio|illumina)}/{chrom,\d+}.cov_greater_than_20.trio_intersect.bed',
            trio_region_bed = 'data/aj_trio/{region}.m{mother_cov}.f{father_cov}.s{son_cov}/trio_covered_regions/{tech,(pacbio|illumina)}/{chrom,\d+}.cov_greater_than_20.trio_intersect.region_filtered.bed',
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

        shell('''
        {BEDTOOLS} intersect -sorted -g {input.genome_file} -a {input.NA24143_bed} -b {input.NA24149_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.parents_bed}
        {BEDTOOLS} intersect -sorted -g {input.genome_file} -a {output.parents_bed} -b {input.NA24385_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.trio_bed}
        {BEDTOOLS} intersect -sorted -g {input.genome_file} -a {output.trio_bed} -b {region_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.trio_region_bed}
        ''')

rule merge_aj_trio_confident:
    params: job_name = 'merge_aj_trio_confident'
    input:  NA24143_bed = 'data/NA24143.hg38/variants/ground_truth/region_filter.bed',
            NA24149_bed = 'data/NA24149.hg38/variants/ground_truth/region_filter.bed',
            NA24385_bed = 'data/NA24385.hg38/variants/ground_truth/region_filter.bed',
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt'
    output: parents_bed = 'data/aj_trio/confident_parent_intersect.bed',
            trio_bed = 'data/aj_trio/confident_trio_intersect.bed',
    shell:
        '''
        {BEDTOOLS} intersect -a {input.NA24143_bed} -b {input.NA24149_bed} > {output.parents_bed}
        {BEDTOOLS} intersect -a {output.parents_bed} -b {input.NA24385_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.trio_bed}
        '''

rule merge_aj_trio_nonconfident:
    params: job_name = 'merge_aj_trio_nonconfident'
    input:  NA24143_bed = 'data/NA24143.hg38/variants/ground_truth/outside_region_filter.bed',
            NA24149_bed = 'data/NA24149.hg38/variants/ground_truth/outside_region_filter.bed',
            NA24385_bed = 'data/NA24385.hg38/variants/ground_truth/outside_region_filter.bed',
            genome_file = 'genome_tracks/hg38.chrom.sizes.natural_order.txt'
    output: parents_bed = 'data/aj_trio/nonconfident_parent_intersect.bed',
            trio_bed = 'data/aj_trio/nonconfident_trio_intersect.bed',
    shell:
        '''
        {BEDTOOLS} intersect -a {input.NA24143_bed} -b {input.NA24149_bed} > {output.parents_bed}
        {BEDTOOLS} intersect -a {output.parents_bed} -b {input.NA24385_bed} | \
        {BEDTOOLS} sort -faidx {input.genome_file} -i - > {output.trio_bed}
        '''

rule generate_coverage_bed_pacbio:
    params: job_name = 'generate_coverage_bed_pacbio.cov_greater_than_20.{id}.{build}.{chrom}.{cov}x'
    input:  'data/{id}.{build}/aligned_reads/pacbio/pacbio.blasr.{chrom}.{cov}x.bam'
    output: 'data/{id}.{build}/aligned_reads/pacbio/pacbio.{chrom,(\d+)}.{cov}x.cov_greater_than_20.bed',
    run:
        maxcov = int(wildcards.cov)*2 # exclude regions more than twice the mean coverage
        shell('''
        {SAMTOOLS} view -F 3844 -q 30 {input} chr{wildcards.chrom} -hb | \
        {BEDTOOLS} genomecov -bga -ibam - | \
        awk '$4 > 20 && $4 < {maxcov}' | \
        {BEDTOOLS} sort -i - | \
        {BEDTOOLS} merge -i - | \
        {BEDTOOLS} sort -i - | \
        ''')

rule generate_coverage_bed_illumina:
    params: job_name = 'generate_coverage_bed_illumina.cov_greater_than_20.{id}.{build}.{chrom}.{cov}x'
    input:  'data/{id}.{build}/aligned_reads/illumina/illumina.{cov}x.bam'
    output: 'data/{id}.{build}/aligned_reads/illumina/illumina.{chrom,(\d+)}.{cov}x.cov_greater_than_20.bed',
    run:
        maxcov = int(wildcards.cov)*2 # exclude regions more than twice the mean coverage
        shell('''
        {SAMTOOLS} view -F 3844 -q 30 {input} chr{wildcards.chrom} -hb | \
        {BEDTOOLS} genomecov -bga -ibam - | \
        awk '$4 > 20 && $4 < {maxcov}' | \
        {BEDTOOLS} sort -i - | \
        {BEDTOOLS} merge -i - | \
        {BEDTOOLS} sort -i - | \
        ''')
