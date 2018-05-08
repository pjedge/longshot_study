import sys
import gzip

def replace_empty_gt_with_reference(input_vcfgz, output_vcf):

    with gzip.open(input_vcfgz, 'rt') as inf, open(output_vcf,'w') as outf:
        for line in inf:
            line = line.strip()

            if line[0] == '#':
                print(line, file=outf)
            else:
                el = line.strip().split('\t')
                if len(el) < 10:
                    continue

                assert(el[8] == 'GT:PS:GQ')
                for i in range(9,len(el)):
                    if el[i] == '.':
                        el[i] = '0/0:.:100.00'

                print('\t'.join(el), file=outf)

#replace_empty_gt_with_reference('data/AK_trio/duplicated_regions/trio_shared_variant_sites/merged/1.with_empties.vcf.gz', 'data/AK_trio/duplicated_regions/trio_shared_variant_sites/merged/1.vcf.gz')
