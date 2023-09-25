import sys
from collections import defaultdict
import re
import bisect

def cigar2pos(start_cigar): 
	pos, cigar = start_cigar.split("_") # start_cigar #
	pos = int(pos)
	if cigar=="93M":
		return([pos], [pos+92]) # 93
	starts = list(); ends = list() 
	cigar_ops = pattern.findall(cigar) 
	if cigar_ops[0][1] == 'S':
		cigar_ops = cigar_ops[1:]
	for length, op in cigar_ops:
		length = int(length)
		if op == 'M':
			starts.append(pos) 
			ends.append(pos+length-1) #-1
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
    for starts_ends in startss_endss:
        for start in starts_ends[0]:
            bisect.insort(starts, start)
        for end in starts_ends[1]:
            bisect.insort(ends, end)
    return find_intersect(starts, ends, th)

def proc(chrom):
	global out 
	for bc in out: 
		for umi in list(out[bc].keys()): # changed size bc of +out[bc][chrom]
			if len(out[bc][umi]) > 1:
				bcumis[bc][umi] = False 
				continue
			strand = list(out[bc][umi].keys())[0]
			if out[bc][umi][strand][0] < dup_min:
				bcumis[bc][umi] = False
				continue
			starts_ends_dups = list()
			start_cigars = out[bc][umi][strand][1]
			start_cigars = start_cigars.split(',')#[1::] ONLY IF out defaultdict str(), ','.joni(out, first start_cigar) => ",.." => empty first, but we switched to normal dict!
			for start_cigar in start_cigars:
				starts_ends_dup = cigar2pos(start_cigar)
				starts_ends_dups.append(starts_ends_dup)
			starts_ends_umi = íntersect(starts_ends_dups, th = dup_min ) #-1
			if starts_ends_umi:
				if chrom in out[bc]:
					out[bc][chrom][umi] = (strand, starts_ends_umi)
				else:
					out[bc][chrom] = {umi:(strand, starts_ends_umi)}
			else:
				bcumis[bc][umi] = False

		if chrom not in out[bc]:
			continue

		if len(out[bc][chrom]) < cov_min :
			for umi in  out[bc][chrom]:
				bcumis[bc][umi] = False
			del out[bc][chrom]
			continue

		starts_ends_umis = defaultdict(lambda: list())
		for umi in  out[bc][chrom]:
			strand, starts_ends_umi = out[bc][chrom][umi]
			starts_ends_umis[strand].append(starts_ends_umi)
		starts_ends_chrbc = dict()
		for strand in starts_ends_umis:
			starts_ends_chrbc[strand] = íntersect(starts_ends_umis[strand], th=cov_min) #-1
		for umi in list(out[bc][chrom].keys()): 
			strand, starts_ends_umi = out[bc][chrom][umi]
			if not starts_ends_chrbc[strand]:
				bcumis[bc][umi] = False
				del out[bc][chrom][umi]
				continue
			coords = íntersect([starts_ends_chrbc[strand], starts_ends_umi], th=2)
			if coords:
				out[bc][chrom][umi] = (strand, coords)
			else:
				bcumis[bc][umi] = False
				del out[bc][chrom][umi]

		if len(out[bc][chrom])<cov_min:
			for umi in  out[bc][chrom]:
				bcumis[bc][umi] = False
			del out[bc][chrom]
			continue

		for umi in  out[bc][chrom]:
			strand, coords = out[bc][chrom][umi]
			bcumis[bc][umi] = ' '.join([chrom, strand, str(coords[0]), str(coords[1])])

		del out[bc][chrom]

def add2out(bc, umi, strand, start_cigar):
	global out
	if bc in out:
		if umi in out[bc]:
			if strand in out[bc][umi]:
				out[bc][umi][strand][0] += 1
				out[bc][umi][strand][1] = ','.join([out[bc][umi][strand][1], start_cigar])
			else:
				bcumis[bc][umi] = False
		else:
			out[bc][umi] = {strand:[1, start_cigar]}
	else:
		out[bc] = {umi:{strand:[1, start_cigar]}}


pattern = re.compile(r'(\d+)([MIDNS])') 

dup_min = int(sys.argv[1]) #5
cov_min = int(sys.argv[2]) #5

filt_bcs=sys.argv[3]
filt_bcs_d=dict()
with open(filt_bcs) as f:
	bc_ind = 1 #0
	for line in f.readlines():
		bc = line.strip()
		filt_bcs_d[bc] = bc_ind


out_file = sys.argv[4]


d_let = {'A':'1','C':'2','G':'3','T':'4'}
d_let_inv = {v:k for k, v in d_let.items()}

bcumis = defaultdict(lambda: defaultdict(lambda: dict() ))
#bcumis = dict()
out = dict()
#out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [int(), str()] ))) #list() # [int(), str()] #, tuple() 

chrom_old=''

for line in sys.stdin:
	line = line.strip().split() #"\t"
	if (line[4]!="255" or int(line[1])>=256 or line[-2]=="CB:Z:-" or line[-1]=="UB:Z:-" or line[2]=="chrM"):
		continue

	bc = line[-2][5::] 
	if bc not in filt_bcs_d:
		continue

	umi = line[-1][5::] 
	if (bc in bcumis) and (umi in bcumis[bc]):
		bcumis[bc][umi] = False
		continue

	strand = line[1][0] ; chrom = line[2]
	start_cigar = '_'.join([line[3], line[5]])

	if chrom == chrom_old:
		add2out(bc, umi, strand, start_cigar)

	elif chrom_old:
		proc(chrom_old)
		out = {bc:{umi:{strand:[1, start_cigar]}}}
		# not add2out() bc we need to clear the variable
		chrom_old = chrom 

	else: # for the first round (chrom_old == False)
		add2out(bc, umi, strand, start_cigar)
		chrom_old = chrom 

proc(chrom)

with open(out_file, 'w') as f:
	for bc in bcumis:
		for umi in bcumis[bc]:
			if bcumis[bc][umi]:
				print(bc, umi, bcumis[bc][umi], file=f)
exit()