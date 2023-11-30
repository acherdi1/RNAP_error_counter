import sys
from collections import defaultdict


subs_file = sys.argv[1]
genes_file = sys.argv[2]

#genes = dict()
default_dict = lambda: defaultdict(lambda: defaultdict())
genes = defaultdict(default_dict)

with open(genes_file) as gf:
	for gline in gf:
		chrom, strand, start, end, enscode, gene = gline.strip().split()
		genes[(chrom, strand)][(int(start), int(end))] = (enscode, gene)

assigned = list()
with open(subs_file) as sf:
	for sline in sf:
		sline = sline.strip()
		_, chrom, strand, pos, _,_ = sline.split()
		pos = int(pos)
		for startend in genes[(chrom, strand)]:
			if (pos >= startend[0] and pos <= startend[1]):
				assigned.append(genes[(chrom, strand)][startend][1])
				#assigned = ','.join([assigned,genes[(chrom, strand)][startend][1]])
		assigned = ','.join(assigned)
		print(sline, assigned)
		assigned = list()
