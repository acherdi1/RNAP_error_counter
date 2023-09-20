#### SAVE COUNTED CIGAR AND TRANSMIT TO THE NEXT STEP!!!

import sys
from collections import defaultdict
import re


def cigar2pos(start, cigar): 
	pos = int(start)
	starts = list() #tuple()
	ends = list() #tuple()
	cigar_ops = pattern.findall(cigar) 
	if cigar_ops[0][1] == 'S':
		cigar_ops = cigar_ops[1:]
	for length, op in cigar_ops:
		length = int(length)
		if op == 'M':
			starts.append(pos) #+= (pos, )
			ends.append(pos+length) #+= (pos + length, )
		elif op == 'I':
			continue
		pos += length
	return(starts, ends) 

def covv(startss, endss, th): 
	startss = sorted(startss, reverse=True)
	endss = sorted(endss, reverse=True)
	starts = startss.pop()
	ends = endss.pop()
	cov = 0
	chunk = False
	cov5 = list() #tuple()
	while startss:
		if starts < ends:
			cov += 1
			if cov >= th and not chunk: #switch_cov5:
				chunk = [starts] #(starts,)
			starts = startss.pop()
		elif starts > ends:
			cov -= 1
			if chunk and cov < th : #switch_cov5
				chunk.append(ends) #+= (ends,)
				cov5.append(chunk) #+= (chunk) #cov5 += chunk
				chunk = False
			ends = endss.pop()
		elif starts == ends:
			starts = startss.pop()
			ends = endss.pop()
		else:
			print("error", starts, ends)
			exit()

	while cov >= th-1:
		if starts < ends:
			cov += 1
			if cov >= th and not chunk: 
				chunk = (starts,)
			starts=endss[0]
		cov -= 1
		if chunk and cov < th : 
			chunk += (ends,)
			cov5 += (chunk, )
			chunk = False
		ends = endss.pop()
	return(cov5)

pattern = re.compile(r'(\d+)([MIDNS])') 


#th = int(sys.argv[1])
dup_min = int(sys.argv[1]) #5
cov_min = int(sys.argv[2]) #5

filt_bcs=sys.argv[2]
filt_bcs_d=dict()
with open(filt_bcs) as f:
	ind = 1
	for line in f.readlines():
		bc = line.strip()
		filt_bcs_d[bc] = ind
		ind += 1


bcumis = dict()
chroms = dict()
chrom_ind = 1

out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))) #, tuple()
chrom_old=''
startss = list() 
endss = list()
covd = dict()

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

			start = int(line[3]); cigar = line[5]

			if chrom == chrom_old:
				out[bc][umi][strand] += 1 
			elif chrom_old:
				for _bc in list(out.keys()):
					for _umi in out[_bc]:
						for _strand in out[_bc][_umi]:
							if out[_bc][_umi][_strand] > dup_min:
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
			if out[_bc][_umi][_strand] > dup_min:
				if (_bc, _umi) not in bcumis:
					bcumis[(_bc, _umi)] = (chrom_old, _strand)
				else:
					bcumis[(_bc, _umi)] = False
