from filter_SNVs import filter_longshot_VCF_for_haplotype_assessment
import sys

rule combine_haplotype_stats:
    params: job_name = 'combine_haplotype_stats.{ext}{individual}.{build}.{method}.{info}',
    input: expand('data/{{individual}}.{{build}}/{{method}}_haplotypes/hap_statistics/{{info}}.{c}{{ext}}',c = chroms)
    output: 'data/{individual}.{build}/{method}_haplotypes/hap_statistics/{info}.all{ext}'
    run:
        err = chs.error_result()
        for f in input:
            err += pickle.load(open(f,'rb'))
        print(err, file=sys.stderr)
        print("switch + mismatch rate: {}".format(err.get_switch_mismatch_rate()), file=sys.stderr)
        pickle.dump(err, open(output[0], "wb" ))

rule haplotype_accuracy_longshot_unfiltered:
    params: job_name = 'haplotype_accuracy_longshot_unfiltered.{individual}.{build}.{tech}{rcov}x.{aligner}.{chrom}',
    input: vcf = 'data/{individual}.{build}/variants/longshot.{tech}.{aligner}.{rcov}x._/{chrom}.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf'
    output: pickle = 'data/{individual}.{build}/longshot_haplotypes/hap_statistics/longshot.{tech}.{aligner}.{rcov}x._.{chrom,(\d+)}.unfiltered.p'
    run:
        err = chs.vcf_vcf_error_rate(input.vcf, input.ground_truth, False)
        print(err, file=sys.stderr)
        pickle.dump(err, open(output.pickle,"wb"))

rule haplotype_accuracy_HapCUT2_unfiltered:
    params: job_name = 'haplotype_accuracy_HapCUT2_unfiltered.{individual}.{build}.illumina30.{tech}{rcov}x.{aligner}.{chrom}',
    input: hap = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.{tech}.{aligner}.{rcov}x/haps/{chrom}',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf'
    output: pickle = 'data/{individual}.{build}/HapCUT2_haplotypes/hap_statistics/illumina.{icov,\d+}x.{tech}.{aligner}.{rcov,\d+}x.{chrom,(\d+)}.unfiltered.p'
    run:
        err = chs.hapblock_vcf_error_rate(input.hap, input.vcf, input.ground_truth, False)
        print(err, file=sys.stderr)
        pickle.dump(err, open(output.pickle, "wb"))

rule haplotype_accuracy_longshot:
    params: job_name = 'haplotype_accuracy_longshot.{individual}.{build}.{tech}{rcov}x.{aligner}.{chrom}',
    input: vcf = 'data/{individual}.{build}/variants/longshot.{tech}.{aligner}.{rcov}x._/{chrom}.PQ30_filtered.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf'
    output: pickle = 'data/{individual}.{build}/longshot_haplotypes/hap_statistics/longshot.{tech}.{aligner}.{rcov}x._.{chrom,(\d+)}.p'
    run:
        err = chs.vcf_vcf_error_rate(input.vcf, input.ground_truth, False)
        print(err, file=sys.stderr)
        pickle.dump(err, open(output.pickle,"wb"))

rule haplotype_accuracy_HapCUT2:
    params: job_name = 'haplotype_accuracy_HapCUT2.{individual}.{build}.illumina30.{tech}{rcov}x.{aligner}.{chrom}',
    input: hap = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.{tech}.{aligner}.{rcov}x/haps/{chrom}.pruned',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf'
    output: pickle = 'data/{individual}.{build}/HapCUT2_haplotypes/hap_statistics/illumina.{icov,\d+}x.{tech}.{aligner}.{rcov,\d+}x.{chrom,(\d+)}.p'
    run:
        err = chs.hapblock_vcf_error_rate(input.hap, input.vcf, input.ground_truth, False)
        print(err, file=sys.stderr)
        pickle.dump(err, open(output.pickle, "wb"))

rule haplotype_accuracy_whatshap:
    params: job_name = 'haplotype_accuracy_whatshap.{individual}.{build}.illumina30.{tech}{rcov}x.{aligner}.{chrom}',
    input: hap = 'data/{individual}.{build}/whatshap_haplotypes/illumina.{icov}x.{tech}.{aligner}.{rcov}x/vcf/{chrom}.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf',
    output: pickle = 'data/{individual}.{build}/whatshap_haplotypes/hap_statistics/illumina.{icov,\d+}x.{tech}.{aligner}.{rcov,\d+}x.{chrom,(\d+)}.p'
    run:
        err = chs.vcf_vcf_error_rate(input.hap, input.ground_truth, False)
        print(err, file=sys.stderr)
        pickle.dump(err, open(output.pickle, "wb"))

NA12878_HG19_PLATINUM_GENOMES_URL = 'https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz'
rule get_haplotyping_ground_truth_NA12878_hg19:
    params: job_name = "get_haplotyping_ground_truth_NA12878_hg19",
    output: pg_vcfgz = 'data/NA12878.1000g/variants/platinum_genomes/with_chr.vcf.gz',
            pg_vcf = 'data/NA12878.1000g/variants/platinum_genomes/with_chr.vcf',
            pg_vcf_nochr = 'data/NA12878.1000g/variants/platinum_genomes/all.vcf',
            pg_chroms = expand('data/NA12878.1000g/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf', chrom=chroms)
    shell:
        '''
        wget {NA12878_HG19_PLATINUM_GENOMES_URL} -O {output.pg_vcfgz}
        gunzip -c {output.pg_vcfgz} > {output.pg_vcf}
        grep -P "^#" {output.pg_vcf} > {output.pg_vcf_nochr}
        grep -Pv "^#" {output.pg_vcf} | awk '{{gsub(/^chr/,""); print}}' >> {output.pg_vcf_nochr}
        for i in {{1..22}}; do
            echo "splitting chr${{i}}"
            grep -P "^chr${{i}}\\t" {output.pg_vcf} | awk '{{gsub(/^chr/,""); print}}' > data/NA12878.1000g/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.${{i}}.vcf
        done
        '''

NA12878_HG38_PLATINUM_GENOMES_URL = 'https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz'
rule get_haplotyping_ground_truth_NA12878_hg38:
    params: job_name = "get_haplotyping_ground_truth_NA12878_hg38",
    output: pg_vcfgz = 'data/NA12878.hg38/variants/platinum_genomes/all.vcf.gz',
            pg_vcf = 'data/NA12878.hg38/variants/platinum_genomes/all.vcf',
            pg_chroms = expand('data/NA12878.hg38/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf', chrom=chroms)
    shell:
        '''
        wget {NA12878_HG38_PLATINUM_GENOMES_URL} -O {output.pg_vcfgz}
        gunzip -c {output.pg_vcfgz} > {output.pg_vcf}
        for i in {{1..22}}; do
            echo "splitting chr${{i}}"
            grep -P "^chr${{i}}\\t" {output.pg_vcf} > data/NA12878.hg38/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.${{i}}.vcf
        done
        '''

rule combine_haplotyping_ground_truth_NA24385:
    params: job_name = "combine_haplotyping_ground_truth_NA24385.{chrom}",
    input:  trio_vcf = 'data/NA24385.hg38/variants/ground_truth/separate_chrom/ground_truth.DECOMPOSED.SNVs_ONLY.{chrom}.vcf',
            tenX_vcf = 'data/NA24385.hg38/variants/ground_truth/separate_chrom/ground_truth.10X.{chrom}.vcf'
    output: vcf = 'data/NA24385.hg38/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf'
    shell: 'python2 merge_phasedvcfs.py {input.trio_vcf} {input.tenX_vcf} > {output.vcf}'

NA24385_10X_VCF_URL = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/NA24385_GRCh38/NA24385_GRCh38_phased_variants.vcf.gz'
rule get_haplotyping_ground_truth_NA24385_10X:
    params: job_name = "get_haplotyping_ground_truth_NA24385_10X",
    output: tenX_vcfgz = 'data/NA24385.hg38/variants/tenX_genomics/all.vcf.gz',
            tenX_vcf = 'data/NA24385.hg38/variants/tenX_genomics/all.vcf',
            tenX_chroms = expand('data/NA24385.hg38/variants/ground_truth/separate_chrom/ground_truth.10X.{chrom}.vcf', chrom=chroms)
    shell:
        '''
        wget {NA24385_10X_VCF_URL} -O {output.tenX_vcfgz}
        gunzip -c {output.tenX_vcfgz} > {output.tenX_vcf}
        for i in {{1..22}}; do
            echo "splitting chr${{i}}"
            grep -P "^chr${{i}}\\t" {output.tenX_vcf} > data/NA24385.hg38/variants/ground_truth/separate_chrom/ground_truth.10X.${{i}}.vcf
        done
        '''

rule prune_longshot_vcf:
    params: job_name = "prune_longshot_haplotype.{individual}.{build}.longshot.{tech}.{aligner}.{rcov}x._.{chrom}",
    input: 'data/{individual}.{build}/variants/longshot.{tech}.{aligner}.{rcov}x._/{chrom}.vcf'
    output: 'data/{individual}.{build}/variants/longshot.{tech}.{aligner}.{rcov}x._/{chrom}.PQ30_filtered.vcf'
    run:
        filter_longshot_VCF_for_haplotype_assessment(input[0], output[0], min_phase_qual=30)

rule prune_HapCUT2_haplotype:
    params: job_name = "prune_HapCUT2_haplotype.{individual}.{build}.illumina{icov}x.{tech}.{aligner}.{rcov}x.{chrom}",
    input: 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.{tech}.{aligner}.{rcov}x/haps/{chrom}'
    output: 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.{tech}.{aligner}.{rcov}x/haps/{chrom}.pruned'
    run:
        ph.prune_hapblock_file(input[0], output[0], snp_conf_cutoff=30.0, split_conf_cutoff=-1.0, use_refhap_heuristic=False)

rule separate_ground_truth_chrom:
    params: job_name = 'separate_ground_truth_chrom.{individual}.{build}.{chrom}',
    input: 'data/{individual}.{build}/variants/ground_truth/ground_truth.DECOMPOSED.SNVs_ONLY.vcf.gz'
    output: 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.DECOMPOSED.SNVs_ONLY.{chrom,(\d+)}.vcf'
    run:
        w_chrom = chr_prefix(wildcards.chrom,wildcards.build)
        shell('''gunzip -c {input} | grep -P '^{w_chrom}\\t' > {output}''')

rule HapCUT2:
    params: job_name = 'HapCUT2.{individual}.{build}.illumina{icov}x.{tech}{rcov}x.{aligner}.{chrom}',
    input: frag = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.{tech}.{aligner}.{rcov}x/fragments/{chrom}',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf'
    output: hap = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov,\d+}x.{tech}.{aligner}.{rcov,\d+}x/haps/{chrom,(\d+)}',
    shell: '{HAPCUT2} --fragments {input.frag} --vcf {input.vcf} --output {output.hap}'

rule extractHAIRS:
    params: job_name = 'extractHAIRS.{individual}.{build}.illumina{icov}.{tech}{rcov}.{aligner}.{chrom}',
    input: bam = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{rcov}x.bam',
           bai = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{rcov}x.bam.bai',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
           hg19 = 'data/genomes/hg19.fa',
           hg19_ix = 'data/genomes/hg19.fa.fai',
           hs37d5 = 'data/genomes/hs37d5.fa',
           hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
           hg38 = 'data/genomes/hg38.fa',
           hg38_ix = 'data/genomes/hg38.fa.fai'
    output: frag = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov,\d+}x.{tech,(pacbio|illumina)}.{aligner}.{rcov,\d+}x/fragments/{chrom,(\d+)}',
    run:
        w_ref = ref_file[wildcards.build]
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        if wildcards.individual == 'NA12878' and wildcards.tech == 'pacbio':
            w_ref = input.hg19 # need to use hg19 rather than 1000g for NA12878...
            chr_prefix_vcf = input.vcf[:-4] + '.chr_prefix.vcf'
            shell('''
                cat {input.vcf} | awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}'> {chr_prefix_vcf};
                {EXTRACTHAIRS} --region chr{w_chrom} --ref {w_ref} --bam {input.bam} --vcf {chr_prefix_vcf} --out {output.frag} --pacbio 1
                ''')
        else:
            tech_flag = ''
            if wildcards.tech == 'pacbio':
                tech_flag = '--ref {w_ref} --pacbio 1'
            elif wildcards.tech == 'ont':
                tech_flag = '--ref {w_ref} --ont 1'

            shell('{EXTRACTHAIRS} ' + tech_flag + ' --region {w_chrom} --bam {input.bam} --vcf {input.vcf} --out {output.frag}')

hg19_chroms = set(['chr{}'.format(i) for i in range(1,23)])
hs37d5_chroms = set([str(i) for i in range(1,23)])
def remove_chr_from_vcf(in_vcf, out_vcf):
    with open(in_vcf, 'r') as inf, open(out_vcf, 'w') as outf:
        for line in inf:
            if line[0] == '#':
                print(line.strip(),file=outf)
                continue
            el = line.strip().split('\t')
            assert(el[0] in hg19_chroms)
            el[0] = el[0][3:]
            assert(el[0] in hs37d5_chroms)
            print("\t".join(el),file=outf)

rule whatshap_phase:
    params: job_name = 'whatshap_phase.{individual}.{build}.illumina{icov}.{tech}{rcov}.{aligner}.{chrom}',
    input: bam = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{rcov}x.bam',
           bai = 'data/{individual}.{build}/aligned_reads/{tech}/{tech}.{aligner}.all.{rcov}x.bam.bai',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
           hg19 = 'data/genomes/hg19.fa',
           hg19_ix = 'data/genomes/hg19.fa.fai',
           hs37d5 = 'data/genomes/hs37d5.fa',
           hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
           hg38 = 'data/genomes/hg38.fa',
           hg38_ix = 'data/genomes/hg38.fa.fai'
    output: vcf = 'data/{individual}.{build}/whatshap_haplotypes/illumina.{icov,\d+}x.{tech,(pacbio|illumina)}.{aligner}.{rcov,\d+}x/vcf/{chrom,(\d+)}.vcf',
    run:
        w_ref = ref_file[wildcards.build]
        w_chrom = chr_prefix(wildcards.chrom, wildcards.build)
        if wildcards.individual == 'NA12878' and wildcards.tech == 'pacbio':
            w_ref = input.hg19 # need to use hg19 rather than 1000g for NA12878...
            chr_prefix_vcf = input.vcf[:-4] + '.chr_prefix.vcf'
            shell('''
                cat {input.vcf} | awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}'> {chr_prefix_vcf};
                whatshap phase --ignore-read-groups --chromosome chr{w_chrom} --reference {input.hg19} -o {output.vcf}.tmp {chr_prefix_vcf} {input.bam}
                ''')
            remove_chr_from_vcf(output.vcf+'.tmp',output.vcf)
        else:
            shell('whatshap phase --ignore-read-groups --chromosome {w_chrom} --reference {w_ref} -o {output.vcf} {input.vcf} {input.bam}')

rule haplotyping_filters:
    params: job_name = 'haplotyping_filters.{individual}.{build}',
    input:  vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/{chrom}.vcf.gz',
            tbi   = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/{chrom}.vcf.gz.tbi'
    output: vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf.gz',
            vcf   = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
    shell:
        '''
        {RTGTOOLS} RTG_MEM=12g vcffilter -k "PASS","." -g 50 -i {input.vcfgz} -o {output.vcfgz};
        gunzip -c {output.vcfgz} > {output.vcf}
        '''

#rule haplotyping_filters:
#    params: job_name = 'haplotyping_filters.{individual}.{build}',
#    input:  vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/{chrom}.vcf.gz',
#            tbi   = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/{chrom}.vcf.gz.tbi'
#    output: tmp_vcfgz = temp('data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/{chrom}.no_quadalleles.vcf.gz'),
#            tmp_tbi   = temp('data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.filtered/{chrom}.no_quadalleles.vcf.gz.tbi'),
#            vcfgz = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf.gz',
#            vcf   = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
#    shell:
#        '''
#        gunzip -c {input.vcfgz} | python3 filter_quadrallelic_sites.py | bgzip -c > {output.tmp_vcfgz}
#        tabix -p vcf {output.tmp_vcfgz}
#        {RTGTOOLS} RTG_MEM=12g vcffilter -k "PASS","." -g 50 -i {output.tmp_vcfgz} -o {output.vcfgz};
#        gunzip -c {output.vcfgz} > {output.vcf}
#        '''
