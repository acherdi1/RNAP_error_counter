
from collections import defaultdict
import sys
import re
#import numpy as np

#default_genome = lambda: defaultdict(lambda: defaultdict(
#	lambda: defaultdict())) # lambda: {1: 0, 2: 0, 3: 0, 4: 0, 0: 0}
#out = defaultdict(default_genome)
#for i in range(seq_len):
#		out[bc][chrom][umi][coords[i]][seq[i]] += 1

#out = defaultdict(lambda: dict())
out = dict()

chroms = dict()
chrom_ind = 0

with open('interm4') as f:
	for line in f: 
		#print(line.strip())
		line = line.strip()
		info, starts, ends = line.split("[")
		bc, chrom, strand, umi = info.split()
		strand = int(strand[0])
		if chrom not in chroms:
			chroms[chrom] = chrom_ind
			chrom_ind += 1
		chrom = chroms[chrom]
		#if chrom != "chr19":
		#	continue
		starts = [int(i) for i in starts.strip('] ').split(', ')]
		ends = [int(i) for i in ends.strip('] ').split(', ')]

		#if (bc, chrom,strand) not in out:
		#	out[bc,chrom,strand] = dict()
		#if umi not in out[bc, chrom,strand]:
		#	out[bc,chrom,strand][umi] = dict()
		#for chunk in range(len(starts)):
		#	for pos in range(starts[chunk], ends[chunk]):
		#		out[bc,chrom,strand][umi][pos] = np.zeros(4,np.uint8) #{1: 0, 2: 0, 3: 0, 4: 0, 0: 0} #[nt[i]] = 1
		if (bc, chrom,strand) not in out:
			out[bc,chrom,strand] = {umi:(starts, ends)}
		else:
			out[bc,chrom,strand][umi] = (starts, ends)
		#print(out)
		#print(out[bc])
		#print(out[bc][chrom])
		#print(out[bc][chrom][umi])
		#print(out[bc][chrom][umi][pos])
		#exit()
		#bcumis.append()

#exit()


chrom_old=''
#for i in chroms:
#	print(i, chroms[i])
#exit()
pattern = re.compile(r'(\d+)([MIDNS])') 
nt2int = {"N":0, "A":1,"C":2,"T":3,"G":4}

for line in sys.stdin:
	#print(line)
	#exit()
	line = line.strip().split() #"\t"
	if line[-2]=="CB:Z:-" or line[-1]=="UB:Z:-" or line[4]!="255" or line[2]=="chrM" or int(line[1])>=256:
		continue
	chrom = line[2]
	if chrom not in chroms:
		continue
	chrom = chroms[chrom]
	bc = line[-2][5::] ; umi = line[-1][5::] ; strand = int(line[1][0])

	if ((bc, chrom,strand) not in out) or (umi not in out[bc,chrom,strand]):
		continue

	#out[bc,chrom,strand][umi] = (starts, ends)
	if isinstance(out[bc,chrom,strand][umi], tuple):
		starts, ends = out[bc,chrom,strand][umi]
		out[bc,chrom,strand][umi] = dict()
		#exit()
		for chunk in range(len(starts)):
			for pos in range(starts[chunk], ends[chunk]):
				out[bc,chrom,strand][umi][pos] = [0,0,0,0,0] #{1: 0, 2: 0, 3: 0, 4: 0, 0: 0} #[nt[i]] = 1

	start = int(line[3]); cigar = line[5]; seq = list(line[9])
	#exit()
	pos = int(start)
	pos_in = 0
	cigar_ops = pattern.findall(cigar) 
	if cigar_ops[0][1] == 'S':
		pos_in += int(cigar_ops[0][0])
		cigar_ops = cigar_ops[1:]
	for length, op in cigar_ops:
		length = int(length)
		if op == 'M':
			#coords[pos_in:(pos_in+length)] = range(pos, pos + length)
			for i in range(length):
				if pos+i in out[bc,chrom,strand][umi]:
					nt_ind = nt2int[seq[pos_in+i]]
					out[bc,chrom,strand][umi][pos+i][nt_ind] += 1
			pos_in += length
		elif op == 'I':
			pos_in += length
			continue
		pos += length
	#print(bc,chrom,strand, umi)
	#print(out[bc,chrom,strand][umi])
	#exit()

exit()



# AGGCCACAGGTGCTAG 0 1 TGCGATGTTAAG                                        
# {28277446: [0, 0, 0, 0, 0], 28277447: [0, 0, 0, 0, 0], 
# 28277448: [0, 0, 0, 0, 0], 28277449: [0, 0, 0, 0, 0], 
# 28277450: [0, 0, 0, 0, 0], 28277451: [0, 0, 0, 0, 0], 
# 28277452: [0, 0, 0, 0, 0], 28277453: [0, 0, 0, 0, 0]}
# 28277345        91M2S
# 28277361        93M
# 28277422        63M1D30M
# 28277431        93M
# 28277432        93M
# 28277507        93M