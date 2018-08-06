
import re
from math import log10

def filter_SNVs(infile, outfile, max_dp, density_count=10, density_len=500, density_qual=50):
    dp_pat = re.compile("DP=(\d+)")

    with open(outfile,'w') as outf:
        lines = []

        with open(infile,'r') as inf:
            for line in inf:
                if line[0] == '#':
                    print(line.strip(),file=outf)
                else:
                    if len(line) < 3:
                        continue
                    el = line.strip().split('\t')
                    chrom = el[0]
                    pos = int(el[1])
                    #qual = float(el[5])
                    format = el[8].split(':')
                    sample = el[9].split(':')
                    gq = None
                    for (tag,data) in zip(format,sample):
                        if tag == 'GQ':
                            gq = float(data)
                    assert(gq != None)
                    lines.append(((chrom,pos,gq),line.strip()))

        filt = [0] * len(lines)
        dp_count = 0
        for i in range(len(lines)):

            depth = int(float(re.findall(dp_pat,lines[i][1])[0]))
            if depth > max_dp:
                filt[i] = 1
                dp_count += 1

            j = i+1
            d = 0
            if lines[i][0][2] > density_qual:
                d += 1
            while j < len(lines):
                if lines[i][0][0] != lines[j][0][0]:
                    break

                if lines[j][0][1] - lines[i][0][1] >= density_len:
                    break

                if lines[j][0][2] > density_qual:
                    d += 1


                if d >= density_count:
                    for k in range(i,j+1):
                        filt[k] = 1

                j += 1

        print("{} variants filtered due to depth".format(dp_count))
        print("{} variants filtered due to density".format(sum(filt)-dp_count))
        #filtered_lines = [l for ((chrom,pos,qual),l),f in zip(lines,filt) if f == 0]

        #for line in filtered_lines:
        #    print(line,file=outf)
        for ((chrom,pos,qual),l),f in zip(lines,filt):
            if f: # filtered out
                el = l.strip().split('\t')
                el[6] = 'fail'
                line = '\t'.join(el)
                print(line,file=outf)
            else:
                print(l,file=outf)


def addlogs(a,b):
    if a > b:
        return a + log10(1.0 + pow(10.0, b - a))
    else:
        return b + log10(1.0 + pow(10.0, a - b))

def filter_reaper_VCF_for_haplotype_assessment(infile, outfile, min_phase_qual=30):
    ph_pat = re.compile("PH=(\d+\.\d+),(\d+\.\d+),(\d+\.\d+),(\d+\.\d+);")

    with open(outfile,'w') as outf:

        with open(infile,'r') as inf:
            for line in inf:
                if line[0] == '#':
                    print(line.strip(),file=outf)
                    continue

                if len(line) < 3:
                    continue
                el = line.strip().split('\t')

                # remove variants not meeting filters
                if el[6] != 'PASS':
                    print(line.strip())
                    continue

                ph_search = re.search(ph_pat,line)
                assert(ph_search)
                # log values of the ph values
                ph1 = float(ph_search.group(1))/-10.0
                ph2 = float(ph_search.group(2))/-10.0
                ph3 = float(ph_search.group(3))/-10.0
                ph4 = float(ph_search.group(4))/-10.0

                # take the 3 smallest log qual values
                # these correspond to the least likely phased genotypes
                vals = sorted([ph1,ph2,ph3,ph4])[:-1]

                # the log sum of the non-max phased genotype probabilities
                qual = addlogs(addlogs(vals[0],vals[1]), vals[2])

                # convert the total to phred score
                phred_qual = qual * -10.0

                # this line has too low phase quality, print filtered line to stdout and move on
                if phred_qual < min_phase_qual:
                    print(line.strip())
                    continue

                # this line has high enough phase quality
                print(line.strip(),file=outf)
