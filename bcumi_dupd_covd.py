## CODE IS TO BE REVISED IF THERE ARE > FLAGS(<256) THAN 0 AND 16 

## output bc umi chr? strand? positions which went through filters, or better chunks
#### SAVE COUNTED CIGAR AND TRANSMIT TO THE NEXT STEP!!!

import sys
from collections import defaultdict
import re
import bisect

def cigar2pos(start_cigar): 
	pos, cigar = start_cigar
	starts = list(); ends = list() 
	cigar_ops = pattern.findall(cigar) 
	if cigar_ops[0][1] == 'S':
		cigar_ops = cigar_ops[1:]
	for length, op in cigar_ops:
		length = int(length)
		if op == 'M':
			starts.append(pos) 
			ends.append(pos+length) 
		elif op == 'I':
			continue
		pos += length
	return(starts, ends) 

def find_intersect(starts, ends, th):
    n = len(starts)
    result_starts = []; result_ends = []
    i = 0; j = 0; overlap_count = 0

    while i < n and j < n:
        if starts[i] <= ends[j]:
            overlap_count += 1
            if overlap_count == th:
                result_starts.append(starts[i])
            i += 1
        else:
            overlap_count -= 1
            if overlap_count == th-1:
                result_ends.append(ends[j])
            j += 1

    if result_starts:
        return result_starts, result_ends
    else:
        return

def íntersect(startss_endss, th):
    starts = list(); ends = list()
    for starts_ends in startss_endss:
        for start in starts_ends[0]:
            bisect.insort(starts, start)
        for end in starts_ends[1]:
            bisect.insort(ends, end)
    return find_intersect(starts, ends, th)

def proc():
	global out, bcumis
	for _bc in out:
		for _umi in out[_bc]:
			if len(out[_bc][_umi]) > 1:
				bcumis[_bc][_umi] = False # if 5+ 5- and 5 next chr
				continue
			_strand = list(out[_bc][_umi].keys())[0]
			if len(out[_bc][_umi][_strand]) < dup_min:
				bcumis[_bc][_umi] = False
				continue
			if (_bc not in bcumis) or (_umi not in bcumis[_bc]):
				starts_ends_dups = list()
				for start_cigar in out[_bc][_umi][_strand]:
					starts_ends_dup = cigar2pos(start_cigar)
					starts_ends_dups.append(starts_ends_dup)
				starts_ends_umi = íntersect(starts_ends_dups, th = dup_min-1 )
				if starts_ends_umi:
					bcumis[_bc][_umi] = (_strand, starts_ends_umi)
			else:
				bcumis[_bc][_umi] = False

		if len([i for i in bcumis[_bc].values() if i])<cov_min:
			#del bcumis[_bc]
			for _umi in  bcumis[_bc]:
				bcumis[_bc][_umi] = False
			continue

		starts_ends_umis = defaultdict(lambda: list())
		for _umi in  bcumis[_bc]:
			if not bcumis[_bc][_umi]:
				continue
			_strand, starts_ends_umi = bcumis[_bc][_umi]
			starts_ends_umis[_strand].append(starts_ends_umi)

		starts_ends_chrbc = dict()
		for _strand in starts_ends_umis:
			starts_ends_chrbc[_strand] = íntersect(starts_ends_umis[_strand], th=cov_min)

		for _umi in bcumis[_bc]: 
			if not bcumis[_bc][_umi]:
				continue
			_strand, starts_ends_umi = bcumis[_bc][_umi]
			coords = íntersect(starts_ends_chrbc[_strand], starts_ends_umi, th=1)
			if coords:
				bcumis[_bc][_umi] = (chrom_old, _strand, coords)
			else:
				bcumis[_bc][_umi] = False

		if len([i for i in bcumis[_bc].values() if i])<cov_min:
			for _umi in  bcumis[_bc]:
				bcumis[_bc][_umi] = False
			continue

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

chroms = dict()
chrom_ind = 1

bcumis = defaultdict(lambda: dict() )
out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: list()))) #, tuple()
chrom_old=''

for line in sys.stdin:
	line = line.strip().split() #"\t"
	if not (line[-2]!="CB:Z:-" and line[-1]!="UB:Z:-" and line[4]=="255" and line[2]!="chrM" and int(line[1])<256):
		continue
	bc = line[-2][5::] 
	if bc not in filt_bcs_d:
		continue

	bc = filt_bcs_d[bc] # replacing bc with its index (string to int)
	umi = line[-1][5::]
	strand = int(line[1]) # 0 forward, 16 reverse; could there be other flags with count>th?
	
	chrom = line[2] # below replacing chrom with its index
	if chrom not in chroms:
		chroms[chrom] = chrom_ind
		chrom_ind += 1
	chrom = chroms[chrom]

	start = int(line[3]); cigar = line[5]

	if chrom == chrom_old:
		out[bc][umi][strand].append( (start, cigar) )
	elif chrom_old:
		proc()
		out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: list())))
		out[bc][umi][strand].append( (start, cigar) ) 
		chrom_old = chrom 
	else: # for the first round (chrom_old == False)
		out[bc][umi][strand].append( (start, cigar) ) 
		chrom_old = chrom 

proc()


for bc in bcumis:
	for umi in bcumis[bc]:
		if bcumis[bc][umi];
			chrom, strand, coords = bcumis[bc][umi]
			print(bc, chrom, strand, umi, *coords )