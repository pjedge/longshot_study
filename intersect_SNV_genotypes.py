
def intersect_SNV_genotypes(clairvoyante_vcf, longshot_vcf, output_vcf):

    snvs = set()

    with open(clairvoyante_vcf,'r') as infile:
        for line in infile:
            if line[0] == '#':
                continue

            el = line.strip().split('\t')

            if el[3] not in 'ACGT' or el[4] not in 'ACGT':
                continue

            assert(el[8][:3] == 'GT:')
            gt = (int(el[9][0]),int(el[9][2]))

            if gt == (1,0):
                gt = (0,1)

            if gt not in [(0,1),(1,1)]:
                continue

            alleles = [el[3],el[4]]

            snvs.add((el[0],int(el[1]),alleles[gt[0]],alleles[gt[1]]))

    with open(longshot_vcf,'r') as infile, open(output_vcf,'w') as outfile:
        for line in infile:
            if line[0] == '#':
                print(line,end='',file=outfile)
                continue

            el = line.strip().split('\t')

            if el[3] not in 'ACGT' or el[4] not in 'ACGT':
                continue

            assert(el[8][:3] == 'GT:')
            gt = (int(el[9][0]),int(el[9][2]))

            if gt == (1,0):
                gt = (0,1)

            if gt not in [(0,1),(1,1)]:
                continue

            alleles = [el[3],el[4]]


            if (el[0],int(el[1]),alleles[gt[0]],alleles[gt[1]]) in snvs:
                print(line,end='',file=outfile)
