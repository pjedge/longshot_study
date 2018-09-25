from filter_SNVs import filter_longshot_VCF_for_haplotype_assessment

rule combine_haplotype_stats:
    params: job_name = 'combine_haplotype_stats.{individual}.{build}.{method}.{info}',
    input: expand('data/{{individual}}.{{build}}/{{method}}_haplotypes/hap_statistics/{{info}}.{c}.p',c = chroms)
    output: 'data/{individual}.{build}/{method}_haplotypes/hap_statistics/{info}.all.p'
    run:
        err = chs.error_result()
        for f in input:
            err += pickle.load(open(f,'rb'))
        pickle.dump(err, open(output[0], "wb" ))

rule haplotype_accuracy_longshot:
    params: job_name = 'haplotype_accuracy_longshot.{individual}.{build}.pacbio{pcov}x.{aligner}.{chrom}',
    input: vcf = 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/{chrom}.filtered.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf',
           contig_size_file = 'genome_tracks/{build}.chrom.sizes.txt'
    output: pickle = 'data/{individual}.{build}/longshot_haplotypes/hap_statistics/longshot.pacbio.{aligner}.{pcov}x._.{chrom,(\d+)}.p'
    run:
        err = chs.vcf_vcf_error_rate(input.vcf, input.ground_truth, input.contig_size_file, False)
        pickle.dump(err, open(output.pickle,"wb"))

rule haplotype_accuracy_HapCUT2:
    params: job_name = 'haplotype_accuracy_HapCUT2.{individual}.{build}.illumina30.pacbio{pcov}x.{aligner}.{chrom}',
    input: hap = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.pacbio.{aligner}.{pcov}x/haps/{chrom}.pruned',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf',
           contig_size_file = 'genome_tracks/{build}.chrom.sizes.txt'
    output: pickle = 'data/{individual}.{build}/HapCUT2_haplotypes/hap_statistics/illumina.{icov,\d+}x.pacbio.{aligner}.{pcov,\d+}x.{chrom,(\d+)}.p'
    run:
        err = chs.hapblock_vcf_error_rate(input.hap, input.vcf, input.ground_truth, input.contig_size_file, False)
        pickle.dump(err, open(output.pickle, "wb"))

NA12878_PLATINUM_GENOMES_URL = 'https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz'
rule get_haplotyping_ground_truth_NA12878:
    params: job_name = "get_haplotyping_ground_truth_NA12878",
    output: pg_vcfgz = 'data/NA12878.1000g/variants/platinum_genomes/with_chr.vcf.gz',
            pg_vcf = 'data/NA12878.1000g/variants/platinum_genomes/with_chr.vcf',
            pg_vcf_nochr = 'data/NA12878.1000g/variants/platinum_genomes/all.vcf',
            pg_chroms = expand('data/NA12878.1000g/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf', chrom=chroms)
    shell:
        '''
        wget {NA12878_PLATINUM_GENOMES_URL} -O {output.pg_vcfgz}
        gunzip -c {output.pg_vcfgz} > {output.pg_vcf}
        grep -P "^#" {output.pg_vcf} > {output.pg_vcf_nochr}
        grep -Pv "^#" {output.pg_vcf} | awk '{{gsub(/^chr/,""); print}}' >> {output.pg_vcf_nochr}
        for i in {{1..22}}; do
            echo "splitting chr${{i}}"
            grep -P "^chr${{i}}\\t" {output.pg_vcf} | awk '{{gsub(/^chr/,""); print}}' > data/NA12878.1000g/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.${{i}}.vcf
        done
        '''

rule get_haplotyping_ground_truth_NA24385:
    params: job_name = "get_haplotyping_ground_truth_NA24385.{chrom}",
    input: 'data/NA24385.hg38/variants/ground_truth/separate_chrom/ground_truth.DECOMPOSED.SNVs_ONLY.{chrom}.vcf'
    output:  'data/NA24385.hg38/variants/ground_truth/separate_chrom/ground_truth.for_haplotyping.{chrom}.vcf'
    shell: 'cp {input} {output}'

rule prune_longshot_vcf:
    params: job_name = "prune_longshot_haplotype.{individual}.{build}.longshot.pacbio.{aligner}.{pcov}x._.{chrom}",
    input: 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/{chrom}.vcf'
    output: 'data/{individual}.{build}/variants/longshot.pacbio.{aligner}.{pcov}x._/{chrom}.filtered.vcf'
    run:
        filter_longshot_VCF_for_haplotype_assessment(input[0], output[0], min_phase_qual=30)

rule prune_HapCUT2_haplotype:
    params: job_name = "prune_HapCUT2_haplotype.{individual}.{build}.illumina{icov}x.pacbio.{aligner}.{pcov}x.{chrom}",
    input: 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.pacbio.{aligner}.{pcov}x/haps/{chrom}'
    output: 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.pacbio.{aligner}.{pcov}x/haps/{chrom}.pruned'
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
    params: job_name = 'HapCUT2.{individual}.{build}.illumina{icov}x.pacbio{pcov}x.{aligner}.{chrom}',
    input: frag = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.pacbio.{aligner}.{pcov}x/fragments/{chrom}',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf'
    output: hap = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov,\d+}x.pacbio.{aligner}.{pcov,\d+}x/haps/{chrom,(\d+)}',
    shell: '{HAPCUT2} --fragments {input.frag} --vcf {input.vcf} --output {output.hap}'

rule extractHAIRS:
    params: job_name = 'extractHAIRS.{individual}.{build}.illumina{icov}.pacbio{pcov}.{aligner}.{chrom}',
    input: bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{pcov}x.split_chroms/{chrom}.bam',
           bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.all.{pcov}x.split_chroms/{chrom}.bam.bai',
           vcf = 'data/{individual}.{build}/variants/freebayes.illumina.aligned.{icov}x.haplotyping_filters/{chrom}.vcf',
           hg19 = 'data/genomes/hg19.fa',
           hg19_ix = 'data/genomes/hg19.fa.fai',
           hs37d5 = 'data/genomes/hs37d5.fa',
           hs37d5_ix = 'data/genomes/hs37d5.fa.fai',
           hg38 = 'data/genomes/hg38.fa',
           hg38_ix = 'data/genomes/hg38.fa.fai'
    output: frag = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov,\d+}x.pacbio.{aligner}.{pcov,\d+}x/fragments/{chrom,(\d+)}',
    run:
        w_ref = ref_file[wildcards.build]
        if wildcards.individual == 'NA12878':
            w_ref = input.hg19 # need to use hg19 rather than 1000g for NA12878...
            chr_prefix_vcf = input.vcf[:-4] + '.chr_prefix.vcf'
            shell('''
                cat {input.vcf} | awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}'> {chr_prefix_vcf};
                {EXTRACTHAIRS} --ref {w_ref} --bam {input.bam} --vcf {chr_prefix_vcf} --out {output.frag} --pacbio 1
                ''')
        else:
            shell('{EXTRACTHAIRS} --ref {w_ref} --bam {input.bam} --vcf {input.vcf} --out {output.frag} --pacbio 1')

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
