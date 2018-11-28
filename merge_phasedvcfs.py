import os,sys

"""
use 10X phased VCF to fix phasing errors in trio-phased VCF files
outputs new filtered VCF where discordant sites are removed

author Vikas bansal, dec 8 2017

usage: python merge_phasedvcfs.py chr1.triophased.vcf chr1.10X.vcf  > chr1.10X_trio.vcf

"""

## trio phased VCF 
variants = {};
File = open(sys.argv[1]); n=0;
for line in File: 
	if line[0] == '#': continue;
	var = line.strip().split('\t'); G = var[9].split(':');
	phase = G[0]; 
	variants[(var[0].strip('chr'),var[1],var[3],var[4])] = [phase,var,-1,0,False]; ## last bit is concordance with 10X
	n +=1
File.close();
print >>sys.stderr, "variants read from ",sys.argv[1],n;

blocks = {}; additionalVAR=0;
## 10X or pacbio phased VCF 
File = open(sys.argv[2]); n=0;
for line in File: 
	if line[0] == '#': continue;
	var = line.strip().split('\t'); 
	if var[6] != "PASS": continue; # only use PASS variants
	G0 = var[8].split(':'); bn=-1;
	for i in xrange(len(G0)): 
		if G0[i] == 'PS': bn = i; ## search for PS tag, blockid
	if bn == -1: continue;
	G = var[9].split(':'); 
	phase = G[0]; block = (var[0].strip('chr'),int(G[bn]))
	#print var
	if phase[1] == '|': 
		if block not in blocks: blocks[block] = [0,0,0]; 
		try: 
			snp = variants[(var[0].strip('chr'),var[1],var[3],var[4])]; 
			snp[2] = block[1]; 
			if snp[0] == '0|1' and phase == '0|1': snp[3] = 1; blocks[block][0] +=1; blocks[block][1] +=1;
			elif snp[0] == '1|0' and phase == '1|0': snp[3] = 1; blocks[block][0] +=1; blocks[block][1] +=1;
			elif snp[0] == '1|0' and phase == '0|1': snp[3] = -1; blocks[block][0] +=1;
			elif snp[0] == '0|1' and phase == '1|0': snp[3] = -1; blocks[block][0] +=1;
			elif snp[0] == '0/1' and (phase == '0|1' or phase == '1|0'): snp[3] = 1000; snp[0] = phase; 
		except KeyError: 
			if phase == '0|1' or phase == '1|0': additionalVAR +=1;
			pass 
File.close();
print >>sys.stderr,"10X file has",additionalVAR,"additional phased heterozygous variants that are not in trio VCF";

## check if each block is in phase 0 or 1 with trio VCF
for block in blocks.iterkeys():
	n= blocks[block][0]; m = blocks[block][1];
	if n >=2: 
		if 10*m > 9*n: blocks[block][2] = 1; 
		elif m < 10*n: blocks[block][2] = -1; 
		print >>sys.stderr, block,blocks[block];

newphased=0;
## identify variants that are discordant in phase between trio and 10X data
for var in variants.iterkeys():
	block = (var[0],variants[var][2]) # chrom,blockid pair
	if block[1] == -1 or blocks[block][2] == 0: continue; 
	p = blocks[block][2]; discor = 0;
	if p == 1 and variants[var][3] == -1: discor =1; 
	if p == -1 and variants[var][3] == 1: discor =1; 
	if variants[var][3] ==1000: 	
		if p == -1 and variants[var][0] == '0|1': variants[var][0] = '1|0';
		elif p == -1 and variants[var][0] == '1|0': variants[var][0] = '0|1'; 
		#print >>sys.stderr, "variant can be phased",var,variants[var]
		newphased +=1;
	if discor ==1: 
		variants[var][4] = True; ## discordant
		#print >>sys.stderr, var[0],var[1],var[2],variants[var][0],blocks[block],variants[var][3]
#print >>sys.stderr, "additional variants phased",newphased;

filtered=0; phased=0;		
## output trio phased VCF filtered version
File = open(sys.argv[1]);
for line in File: 
	if line[0] == '#': 
		print line,
		continue;
	var = line.strip().split('\t'); G = var[9].split(':');
	discor = variants[(var[0].strip('chr'),var[1],var[3],var[4])];
	if discor[3] == 1000: ## new SNP phased using 10X data, unphased in trio VCF
		vcfstring = '\t'.join(discor[1][0:8] + ['GT:NT'] + [discor[0] + ':1']);
		#print >>sys.stderr,"need to update";
		#print vcfstring; ### changed this option to not output 11/23/2018
		pass; 
	elif discor[4] == False: phased +=1; print line,
	elif discor[4] == True: filtered +=1; 
File.close();

print >>sys.stderr,"filtered",filtered,"variants in consensus phased VCF","total",phased;
	
