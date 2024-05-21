import sys
from collections import defaultdict


subs_file = sys.argv[1]
genes_file = sys.argv[2]

#genes = dict()
default_dict = lambda: defaultdict(lambda: defaultdict())
genes = defaultdict(default_dict)

with open(genes_file) as gf:
	for gline in gf:
		#chrom, strand, start, end, enscode, gene = gline.strip().split()
		chrom, _, start, end, enscode, gene = gline.strip().split()
		genes[chrom][(int(start), int(end))] = (enscode, gene)

assigned = list()
poss = list()
#switch={'+':'-','-':'+'}
with open(subs_file) as sf:
	for sline in sf:
		sline = sline.strip().split()
		#_, chrom, strand, pos, _,_ = sline.split()
		#_, chrom, _, pos, _,_,_ = sline.split()
		chrom = sline[1] #.lstrip("chr")
		pos = int(sline[3]) #; strand=switch[strand]
		#for startend in genes[(chrom, strand)]:
		for startend in genes[chrom]:
			if (pos >= startend[0] and pos <= startend[1]):
				#assigned.append(genes[(chrom, strand)][startend][1])
				assigned.append('|'.join([str(genes[chrom][startend][1]), str(startend[1]-startend[0])]))
				poss.append( str( (pos - startend[0]) / startend[1] ))
				#assigned = ','.join([assigned,genes[(chrom, strand)][startend][1]])
		assigned = ','.join(assigned); poss = ','.join(poss)

		print(sline[0], chrom, sline[2], pos, sline[4], sline[5], assigned, poss)
		assigned = list(); poss = list()
