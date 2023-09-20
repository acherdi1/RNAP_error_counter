import sys
from collections import defaultdict
th = int(sys.argv[1])

filt_bcs=sys.argv[2]

filt_bcs_d=dict()
with open(filt_bcs) as f:
	ind = 1
	for line in f.readlines():
		bc = line.strip()
		filt_bcs_d[bc] = ind
		ind += 1

chrom_old=''

out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))) #, tuple()

bcumis = dict()
chroms = dict()
chrom_ind = 1

for line in sys.stdin:
	line = line.strip().split() #"\t"
	if (line[-2]!="CB:Z:-" and line[-1]!="UB:Z:-" and line[4]=="255" and line[2]!="chrM" and int(line[1])<256):
		bc = line[-2][5::] 
		#print(chrom_old, chrom, end=' ')
		if bc in filt_bcs_d:

			bc = filt_bcs_d[bc] # replacing bc with its index (string to int)
			umi = line[-1][5::]
			strand = int(line[1]) # 0 forward, 16 reverse; could there be other flags with count>th?
			chrom = line[2]

			# below replacing chrom with its index
			if chrom not in chroms:
				chroms[chrom] = chrom_ind
				chrom_ind += 1
			chrom = chroms[chrom]

			if chrom == chrom_old:
				out[bc][umi][strand] += 1 #; print('==')
			elif chrom_old:
				for _bc in list(out.keys()):
					for _umi in out[_bc]:
						for _strand in out[_bc][_umi]:
							if out[_bc][_umi][_strand]>=th:
								#print(_bc, _umi)
								# stuff below works for strand check as well
								# as long as there are only 0 and 16 flags
								# 0 forward, 16 reverse
								# could there be other flags with count>th?
								if (_bc, _umi) not in bcumis:
									bcumis[(_bc, _umi)] = (chrom_old, _strand)
								else:
									bcumis[(_bc, _umi)] = False

				out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
				out[bc][umi][strand] += 1
				chrom_old = chrom #; print('!=,+')
			else: # for the first round (chrom_old == False)
				out[bc][umi][strand] += 1
				chrom_old = chrom #; print('!=,-')

for _bc in list(out.keys()):
	for _umi in out[_bc]:
		for _strand in out[_bc][_umi]:
			if out[_bc][_umi][_strand]>=th:
				if (_bc, _umi) not in bcumis:
					bcumis[(_bc, _umi)] = (chrom_old, _strand)
				else:
					bcumis[(_bc, _umi)] = False
