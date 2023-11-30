import sys
from collections import defaultdict
# from natsort import natsorted

input_file = sys.argv[1]
default_dict = lambda: defaultdict(lambda: 0)
d = defaultdict(default_dict)


with open(input_file) as f:
	for line in f:
		#print(line)
		line = line.strip().split()
		if len(line) == 7:
			_, chrom, strand, pos, cons, mm, gene = line
			d[(chrom, strand, pos, gene)][(cons,mm)] += 1
		else:
			_, chrom, strand, pos, cons, mm = line
			d[(chrom, strand, pos)][(cons,mm)] += 1

for pos in d: #natsorted(d.keys()):
	#keys=list(d[pos].keys()) #print(keys)
	if len(d[pos])>1: #==2 and keys[0][1]==keys[1][0]: # True: #len(d[pos])>1: #==2 and keys[0][1]==keys[1][0]:
		for cons_mm in d[pos]:
			print(*pos, *cons_mm, d[pos][cons_mm])
		print()
exit()
