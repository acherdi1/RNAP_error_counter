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
    starts.append(ends[-1]+1)

    result_starts = []; result_ends = []
    i = 0; j = 0; overlap_count = 0

    while i <= n and j < n:
        #print(result_starts, result_ends, i,j,n,overlap_count)
        if starts[i] <= ends[j]:
            overlap_count += 1
            if overlap_count == th:
                result_starts.append(starts[i])
            i += 1
        else:
            overlap_count -= 1
            if overlap_count == th-1:
                result_ends.append(ends[j])
                if i == n:
                    break
            j += 1

    if result_starts:
        return result_starts, result_ends
    else:
        return

def íntersect(startss_endss, th):
    starts = list(); ends = list()
    #print('ss',startss_endss)
    for starts_ends in startss_endss:
        #print('s', starts_ends)
        for start in starts_ends[0]:
            bisect.insort(starts, start)
        for end in starts_ends[1]:
            bisect.insort(ends, end)
    #print('in', startss_endss, starts, ends)
    #if len(starts) != len(ends):
    #	print(startss_endss, starts, ends)
    #	exit()
    #print(find_intersect(starts, ends, th))
    #exit()
    return find_intersect(starts, ends, th)

def proc():
	global out, bcumis, chrom_old
	for _bc in out:
		for _umi in out[_bc]:
			if len(out[_bc][_umi]) > 1:
				bcumis[_bc, _umi] = False # if 5+ 5- and 5 next chr
				continue
			_strand = list(out[_bc][_umi].keys())[0]
			if len(out[_bc][_umi][_strand]) < dup_min:
				bcumis[_bc, _umi] = False
				continue
			if (_bc, _umi) not in bcumis:
				starts_ends_dups = list()
				for start_cigar in out[_bc][_umi][_strand]:
					starts_ends_dup = cigar2pos(start_cigar)
					starts_ends_dups.append(starts_ends_dup)
				#print(starts_ends_dups)
				#exit()
				#print("dups2umi")
				starts_ends_umi = íntersect(starts_ends_dups, th = dup_min-1 )
				#print('out', starts_ends_umi)
				if starts_ends_umi:
					bcumis[_bc][chrom_old][_umi] = (_strand, starts_ends_umi)
			else:
				bcumis[_bc, _umi] = False

		if chrom_old not in bcumis[_bc]:
			continue
		#else:
		#	print('1')
		if len(bcumis[_bc][chrom_old])<cov_min:
			#print('len small 1')
			#del bcumis[_bc]
			for _umi in  bcumis[_bc][chrom_old]:
				bcumis[_bc, _umi] = False
			continue

		#len(bcumis[_bc][_umi])==3
		# that's horrific, 
		#after each chr it would check ALL the umis of each BC present,
		# instead of only working with the output from last chr
		# FIXED !! intermediate bcumis[_bc][chrom_old]

		starts_ends_umis = defaultdict(lambda: list())
		for _umi in  bcumis[_bc][chrom_old]:
			#if not bcumis[_bc][_umi] or len(bcumis[_bc][_umi])==3:
			#	continue
			_strand, starts_ends_umi = bcumis[_bc][chrom_old][_umi]
			#print(starts_ends_umi)
			starts_ends_umis[_strand].append(starts_ends_umi)

		starts_ends_chrbc = dict()
		for _strand in starts_ends_umis:
			#print("umis2chrbc")
			#print(starts_ends_umis[_strand])
			starts_ends_chrbc[_strand] = íntersect(starts_ends_umis[_strand], th=cov_min)
		for _umi in list(bcumis[_bc][chrom_old].keys()): 
			#if not bcumis[_bc][_umi]  or len(bcumis[_bc][_umi])==3:
			#	continue
			_strand, starts_ends_umi = bcumis[_bc][chrom_old][_umi]
			#coords = íntersect(starts_ends_chrbc[_strand], starts_ends_umi, th=1)
			#print(starts_ends_chrbc[_strand], starts_ends_umi)
			#exit()
			if not starts_ends_chrbc[_strand]:
				bcumis[_bc,_umi] = False
				del bcumis[_bc][chrom_old][_umi]
				continue
			coords = íntersect([starts_ends_chrbc[_strand], starts_ends_umi], th=1)
			if coords:
				bcumis[_bc][chrom_old][_umi] = (_strand, coords)
			else:
				bcumis[_bc,_umi] = False
				del bcumis[_bc][chrom_old][_umi]


		if len(bcumis[_bc][chrom_old])<cov_min:
			#print('len small 2')
			for _umi in  bcumis[_bc][chrom_old]:
				bcumis[_bc,_umi] = False
		else:
			#print(2)
			for _umi in  bcumis[_bc][chrom_old]:
				bcumis[_bc,_umi] = True #(chrom_old, _strand, coords)
				_bc_ = l_bcs[_bc]; _chr = list(chroms.keys())[chrom_old]
				_umi_ = ''.join([d_let_inv[i] for i in list(str(_umi))])
				print(_bc_, _chr, _strand, _umi_, *coords )
		del bcumis[_bc][chrom_old]


pattern = re.compile(r'(\d+)([MIDNS])') 

#th = int(sys.argv[1])
dup_min = int(sys.argv[1]) #5
cov_min = int(sys.argv[2]) #5

filt_bcs=sys.argv[3]
filt_bcs_d=dict()
with open(filt_bcs) as f:
	ind = 0
	for line in f.readlines():
		bc = line.strip()
		filt_bcs_d[bc] = ind
		ind += 1
l_bcs = list(filt_bcs_d.keys())

chroms = dict()
chrom_ind = 0
d_let = {'A':'1','C':'2','G':'3','T':'4'}
d_let_inv = {v:k for k, v in d_let.items()}

bcumis = defaultdict(lambda: defaultdict(lambda: dict() ))
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
	umi = int(''.join([d_let[i] for i in list(line[-1][5::])]))
	strand = int(line[1][0]) # 0 forward, 16 reverse; could there be other flags with count>th?
	
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
		#if chrom == "chr3":
		#	exit()
		out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: list())))
		out[bc][umi][strand].append( (start, cigar) ) 
		chrom_old = chrom 
	else: # for the first round (chrom_old == False)
		out[bc][umi][strand].append( (start, cigar) ) 
		chrom_old = chrom 

proc()


exit()

l_chrs = list(chroms.keys())
for bc in bcumis:
	for umi in bcumis[bc]:
		if bcumis[bc][umi]:
			chrom, strand, coords = bcumis[bc][umi]
			_bc = l_bcs[bc]; _chr = l_chrs[chrom]
			_umi = ''.join([d_let_inv[i] for i in list(str(umi))])
			print(_bc, _chr, strand, umi, *coords )