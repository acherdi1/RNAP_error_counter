import sys
from collections import defaultdict
#from time import gmtime, strftime

input_file = sys.argv[1] #'DR_age/res/DR1_old/subs3.wg.txt' #
reads_wg_file = sys.argv[2] #${out}/reads.wg.txt
genes_file = sys.argv[3] #genes.gex_mm10_2020_A
out_gene_cov = sys.argv[4]
out_gene_abs = sys.argv[5]

gene_cov = defaultdict(lambda: 0)

with open(reads_wg_file) as f:
	for line in f:
		line = line.strip().split()
		if len(line)==7:
			gene = line[-1]
			gene_cov[gene] += 1

#print('gene_mm done', strftime("%a, %d %b %Y %H:%M:%S", gmtime()))

gene_mm = {gene:0 for gene in gene_cov} #defaultdict(lambda: defaultdict(0))
gene_absent = defaultdict(lambda: 0)
with open(input_file) as f:
	for line in f:
		line = line.strip().split()
		if len(line)==7:
			gene = line[-1]
			if ',' in gene:
				gene = gene.split(',')
				for g in gene:
					if g in gene_mm:
						gene_mm[g] += 1
					else:
						gene_absent[g] += 1
			else:
				if gene in gene_mm:
					gene_mm[gene] += 1
				else:
					gene_absent[gene] += 1

#print('gene_cov done', strftime("%a, %d %b %Y %H:%M:%S", gmtime()))

gene_len = {gene:0 for gene in gene_cov}
gene_info = dict()
with open(genes_file) as f:
	for line in f:
		line = line.strip().split()
		gene = line[-1]
		gene_len[gene] = int(line[3]) - int(line[2])
		gene_info[gene] = (line[0], line[1])

#print('gene_len done', strftime("%a, %d %b %Y %H:%M:%S", gmtime()))

with open(out_gene_cov, "w") as f:
	for gene in gene_cov:
		print(*gene_info[gene], gene, gene_len[gene], gene_cov[gene], gene_mm[gene], file=f)

with open(out_gene_abs, "w") as f:
	for gene in gene_absent:
		print(*gene_info[gene], gene, gene_len[gene], gene_absent[gene], file=f)
