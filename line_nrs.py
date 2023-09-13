import sys

out = dict()
bcumis_mt_th = sys.argv[1]
with open(bcumis_mt_th) as f:
	for line in f.readlines():
		line = line.strip()
		out[line] = [False, list()] #tuple()]

line_nr = 0
for line in sys.stdin:
	line_nr += 1
	line = line.strip().split() #"\t"
	if (line[2]!="chrM" and line[4]=="255" and 
		line[-2]!="CB:Z:-" and line[-1]!="UB:Z:-"):
		bc_umi = ' '.join([line[-2][5::], line[-1][5::]])
		if bc_umi not in out:
			continue
		chrom= ''.join([line[1][0], line[2]])
		if not out[bc_umi][0]:
			out[bc_umi][0] = chrom
		else:
			if out[bc_umi][0] != chrom:
				out[bc_umi][0] = 1
		out[bc_umi][1].append(line_nr) # += (line_nr,)

for bc_umi in list(out.keys()):
	if out[bc_umi][0] != 1:
		for i in out[bc_umi][1]:
			print(i)

