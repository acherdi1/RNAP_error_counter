import sys
from collections import defaultdict
th = int(sys.argv[1])

filt_bcs=sys.argv[2]

filt_bcs_d=dict()
with open(filt_bcs) as f:
	for line in f.readlines():
		bc = line.strip()
		filt_bcs_d[bc] = 1

chrom_old=''

out = defaultdict(lambda: defaultdict(lambda: 0)) #, tuple()

for line in sys.stdin:
	line = line.strip().split() #"\t"
	if (line[-2]!="CB:Z:-" and line[-1]!="UB:Z:-" and line[4]=="255" and line[2]!="chrM" and int(line[1])<256):
		bc = line[-2][5::] 
		#print(chrom_old, chrom, end=' ')
		if bc in filt_bcs_d:
			umi = line[-1][5::]; chrom = line[2]
			if chrom == chrom_old:
				out[bc][umi] += 1 #; print('==')
			elif chrom_old:
				for _bc in list(out.keys()):
					for _umi in out[_bc]:
						if out[_bc][_umi]>=th:
							print(_bc, _umi)
				out = defaultdict(lambda: defaultdict(lambda: 0))
				out[bc][umi] += 1
				chrom_old = chrom #; print('!=,+')
			else:
				out[bc][umi] += 1
				chrom_old = chrom #; print('!=,-')
for _bc in list(out.keys()):
	for _umi in out[_bc]:
		if out[_bc][_umi]>=th:
			print(_bc, _umi)
