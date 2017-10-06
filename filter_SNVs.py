
import re

filter_SNVs(infile, outfile, max_dp, density_count=10, density_len=500, density_qual=50):
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
                    qual = float(el[5])
                    lines.append(((chrom,pos,qual),line.strip()))

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

        print("{} variants filtered due to depth".format(dp_count)
        print("{} variants filtered due to density".format(sum(filt)-dp_count))
        filtered_lines = [l for ((chrom,pos,qual),l),f in zip(lines,filt) if f == 0]

        for line in filtered_lines:
            print(line,file=outf)
