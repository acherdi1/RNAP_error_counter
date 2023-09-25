import sys
from collections import defaultdict
import re
import bisect


# only add2out(), exit when chrom_old => 42s 1.3GB
# StringIO => 26s, 0.3GB!!!

#sys.getsizeof({k:1 for k in range(5000)})
#147552
#sys.getsizeof({k:1 for k in range(50000)})
#2621536
#sys.getsizeof({1:{k:1 for k in range(25000)}, 2:{k:1 for k in range(25000)}})
#232
# so nestedness is a cure against extramem issues


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
			#for umi_p2 in list(out[bc][umi_p1].keys()):
			if len(out[bc][umi]) > 1:
				add2bcumis(bc, umi, False) 
				continue
			strand = list(out[bc][umi].keys())[0]
			if out[bc][umi][strand][0] < dup_min:
				add2bcumis(bc, umi, False)
				continue
			starts_ends_dups = list()
			start_cigars = out[bc][umi][strand][1]
			start_cigars = start_cigars.split(',') #[1::] #ONLY IF out defaultdict str(), ','.joni(out, first start_cigar) => ",.." => empty first, but we switched to normal dict!
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
				add2bcumis(bc, umi, False)

		if chrom not in out[bc]:
			continue

		if len(out[bc][chrom]) < cov_min :
			for umi in  out[bc][chrom]:
				add2bcumis(bc, umi, False)
			del out[bc][chrom]
			continue

		starts_ends_umis = defaultdict(lambda: list()) # [[],[]] #
		for umi in  out[bc][chrom]:
			strand, starts_ends_umi = out[bc][chrom][umi]
			starts_ends_umis[strand].append(starts_ends_umi)
		starts_ends_chrbc = dict() # [None, None] #
		for strand in starts_ends_umis:
			starts_ends_chrbc[strand] = íntersect(starts_ends_umis[strand], th=cov_min) #-1
		
			if not starts_ends_chrbc[strand]:
				for umi in list(out[bc][chrom].keys()): 
					if out[bc][chrom][umi][0] == strand:
						add2bcumis(bc, umi, False)
						del out[bc][chrom][umi]
						continue
		for umi in list(out[bc][chrom].keys()): 
			strand, starts_ends_umi = out[bc][chrom][umi]
			coords = íntersect([starts_ends_chrbc[strand], starts_ends_umi], th=2)
			if coords:
				out[bc][chrom][umi] = (strand, coords)
			else:
				add2bcumis(bc, umi, False)
				del out[bc][chrom][umi]

		if len(out[bc][chrom])<cov_min:
			for umi in  out[bc][chrom]:
				add2bcumis(bc, umi, False)
			del out[bc][chrom]
			continue

		for umi in  out[bc][chrom]:
			strand, coords = out[bc][chrom][umi]
			#bcumis[bc][umi] = ' '.join([chrom, strand, str(coords[0]), str(coords[1])]) #True #
			value = ' '.join([chrom, strand, str(coords[0]), str(coords[1])]) #True #
			add2bcumis(bc, umi, value)

		del out[bc][chrom]

def add2out(bc, umi, strand, start_cigar):
	global out
	if bc in out:
		if umi in out[bc]:
			#if isinstance(out[bc][umi],int):
			#	out[bc][umi] = {strand:[1, start_cigar]}
			if strand in out[bc][umi]:
				out[bc][umi][strand][0] += 1
				out[bc][umi][strand][1] = ','.join([out[bc][umi][strand][1], start_cigar])
			else:
				add2bcumis(bc, umi, False)
		else:
			out[bc][umi] = {strand:[1, start_cigar]}
	else:
		out[bc] = {umi:{strand:[1, start_cigar]}}

def add2bcumis(bc, umi, value):
	global bcumis
	if bc in bcumis:
		bcumis[bc][umi] = value
	else:
		bcumis[bc] = {umi:value}


pattern = re.compile(r'(\d+)([MIDNS])') 

dup_min = int(sys.argv[1]) #5
cov_min = int(sys.argv[2]) #5

filt_bcs=sys.argv[3]
filt_bcs_d=dict()
with open(filt_bcs) as f:
	bc_ind = 1 
	for line in f.readlines():
		bc = line.strip()
		filt_bcs_d[bc] = bc_ind

out_file = sys.argv[4]

#bcumis = defaultdict(lambda: defaultdict(lambda: dict() ))
bcumis = dict()
out = dict()
#out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [int(), str()] ))) #list() # [int(), str()] #, tuple() 
#out_def=defaultdict(lambda: defaultdict(lambda: [int(), str()] ))
#out = defaultdict(lambda: out_def)

chrom_old=''

for line in sys.stdin:
	line = line.strip().split() #"\t"
	if (line[4]!="255" or int(line[1])>=256 or line[-2][5::] not in filt_bcs_d or line[-2]=="CB:Z:-" or line[-1]=="UB:Z:-" or line[2]=="chrM"):
		continue

	bc = line[-2][5::]; umi = line[-1][5::] 
	if (bc in bcumis) and (umi in bcumis[bc]):
		add2bcumis(bc, umi, False)
		continue

	strand = line[1][0] 
	if (bc in out) and (umi in out[bc]) and (strand not in out[bc][umi]): 
		add2bcumis(bc, umi, False) ; continue

	chrom = line[2]
	start_cigar = '_'.join([line[3], line[5]])

	if chrom != chrom_old:
		if chrom_old:
			proc(chrom_old)
			out = dict()
			#out = defaultdict(lambda: out_def)
		chrom_old = chrom 
	add2out(bc, umi, strand, start_cigar)
	#out[bc][umi][strand][0] += 1
	#out[bc][umi][strand][1] = ','.join([out[bc][umi][strand][1], start_cigar])


proc(chrom)

with open(out_file, 'w') as f:
	for bc in bcumis:
		for umi in bcumis[bc]:
			if bcumis[bc][umi]:
				print(bc, umi, bcumis[bc][umi], file=f)
exit()