from filter_SNVs import filter_reaper_VCF_for_haplotype_assessment

rule combine_haplotype_stats:
    params: job_name = 'combine_haplotype_stats.{individual}.{build}.{method}.{info}',
    input: expand('data/{{individual}}.{{build}}/{{method}}_haplotypes/hap_statistics/{{info}}.{c}.p',c = chroms)
    output: 'data/{individual}.{build}/{method}_haplotypes/hap_statistics/{info}.all.p'
    run:
        err = chs.error_result()
        for f in input:
            err += pickle.load(open(f,'rb'))
        pickle.dump(err, open(output, "wb" ))

rule haplotype_accuracy_reaper:
    params: job_name = 'haplotype_accuracy_reaper.{individual}.{build}.pacbio{pcov}x.PQ{GQ}.{aligner}.{chrom}',
    input: vcf = 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{pcov}x.-z/{chrom}.PQ{GQ}.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.DECOMPOSED.SNVs_ONLY.{chrom}.vcf',
           contig_size_file = 'genome_tracks/{build}.chrom.sizes.txt'
    output: pickle = 'data/{individual}.{build}/reaper_haplotypes/hap_statistics/reaper.pacbio.{aligner}.{pcov}x.-z.{chrom,(\d+)}.p'
    run:
        err = chs.vcf_vcf_error_rate(input.vcf, input.ground_truth, input.contig_size_file, False)
        pickle.dump(err, open(output.pickle,"wb"))

rule haplotype_accuracy_HapCUT2:
    params: job_name = 'haplotype_accuracy_HapCUT2.{individual}.{build}.illumina30.pacbio{pcov}x.{aligner}.{chrom}',
    input: hap = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.pacbio.{aligner}.{pcov}x/haps/{chrom}.pruned',
           vcf = 'data/{individual}.{build}/variants/illumina_{icov}x.filtered/{chrom}.vcf',
           ground_truth = 'data/{individual}.{build}/variants/ground_truth/separate_chrom/ground_truth.DECOMPOSED.SNVs_ONLY.{chrom}.vcf'
    output: pickle = 'data/{individual}.{build}/HapCUT2_haplotypes/hap_statistics/illumina.{icov,\d+}x.pacbio.{aligner}.{pcov,\d+}x.{chrom,(\d+)}.p'
    run:
        err = chs.hapblock_vcf_error_rate(input.hap, input.vcf, input.ground_truth, input.contig_size_file, False)
        pickle.dump(err, open(output.pickle, "wb"))

rule prune_reaper_vcf:
    params: job_name = "prune_reaper_haplotype.{individual}.{build}.reaper.pacbio.{aligner}.{pcov}x.-z.{chrom}",
    input: 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{pcov}x.-z/{chrom}.vcf'
    output: 'data/{individual}.{build}/variants/reaper.pacbio.{aligner}.{pcov}x.-z/{chrom}.PQ{GQ}.vcf'
    run:
        filter_reaper_VCF_for_haplotype_assessment(input, output, min_phase_qual=30)

rule prune_HapCUT2_haplotype:
    params: job_name = "prune_HapCUT2_haplotype.{individual}.{build}.illumina{icov}x.pacbio.{aligner}.{pcov}x.{chrom}",
    input: 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.pacbio.{aligner}.{pcov}x/haps/{chrom}'
    output: 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov}x.pacbio.{aligner}.{pcov}x/haps/{chrom}.pruned'
    run:
        prune_hapblock_file(input, output, snp_conf_cutoff=50.0, split_conf_cutoff=-1.0, use_refhap_heuristic=False)

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
           vcf = 'data/{individual}.{build}/variants/illumina_{icov}x.filtered/{chrom}.vcf'
    output: hap = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov,\d+}x.pacbio.{aligner}.{pcov,\d+}x/haps/{chrom,(\d+)}',
    shell: '{HAPCUT2} --frag {input.frag} --vcf {input.vcf} --out {output.hap}'

rule extractHAIRS:
    params: job_name = 'extractHAIRS.{individual}.{build}.illumina{icov}.pacbio{pcov}.{aligner}.{chrom}',
    input: bam = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.{chrom}.{pcov}x.bam',
           bai = 'data/{individual}.{build}/aligned_reads/pacbio/pacbio.{aligner}.{chrom}.{pcov}x.bam.bai',
           vcf = 'data/{individual}.{build}/variants/illumina_{icov}x.filtered/{chrom}.vcf'
    output: frag = 'data/{individual}.{build}/HapCUT2_haplotypes/illumina.{icov,\d+}x.pacbio.{aligner}.{pcov,\d+}x/fragments/{chrom,(\d+)}',
    run:
        w_ref = ref_file[wildcards.build]
        if wildcards.individual == 'NA12878':
            chr_prefix_vcf = input.vcf[:-4] + '.chr_prefix.vcf'
            shell('''cat {input.vcf} | awk '{{if($0 !~ /^#/) print "chr"$0; else print $0}}'> {chr_prefix_vcf}''')
            shell('{EXTRACTHAIRS} --ref {w_ref} --bam {input.bam} --vcf {chr_prefix_vcf} --out {output.frag} --pacbio 1')
        else:
            shell('{EXTRACTHAIRS} --ref {w_ref} --bam {input.bam} --vcf {input.vcf} --out {output.frag} --pacbio 1')
