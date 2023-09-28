import sys
from collections import defaultdict
import re
import bisect
from microdict import mdict
import io 

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
				bad2bcumis(bc, umi) 
				continue
			strand = list(out[bc][umi].keys())[0]
			if out[bc][umi][strand][0] < dup_min:
				bad2bcumis(bc, umi)
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
				bad2bcumis(bc, umi)

		if chrom not in out[bc]:
			continue

		if len(out[bc][chrom]) < cov_min :
			for umi in  out[bc][chrom]:
				bad2bcumis(bc, umi)
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
						bad2bcumis(bc, umi)
						del out[bc][chrom][umi]
						continue
		for umi in list(out[bc][chrom].keys()): 
			strand, starts_ends_umi = out[bc][chrom][umi]
			coords = íntersect([starts_ends_chrbc[strand], starts_ends_umi], th=2)
			if coords:
				out[bc][chrom][umi] = (strand, coords)
			else:
				bad2bcumis(bc, umi)
				del out[bc][chrom][umi]

		if len(out[bc][chrom])<cov_min:
			for umi in  out[bc][chrom]:
				bad2bcumis(bc, umi)
			del out[bc][chrom]
			continue

		for umi in  out[bc][chrom]:
			strand, coords = out[bc][chrom][umi]
			#bcumis[bc][umi] = ' '.join([chrom, strand, str(coords[0]), str(coords[1])]) #True #
			value = ' '.join([chrom, strand, str(coords[0]), str(coords[1]), ";"]) #True #
			good2bcumis(bc, umi, value)

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
				bad2bcumis(bc, umi)
		else:
			out[bc][umi] = {strand:[1, start_cigar]}
	else:
		out[bc] = {umi:{strand:[1, start_cigar]}}


pattern = re.compile(r'(\d+)([MIDNS])') 

dup_min = int(sys.argv[1]) #5
cov_min = int(sys.argv[2]) #5

filt_bcs=sys.argv[3]
bcumis=dict()
keeper= io.StringIO() #[False]
keeper.write(';')
keeper_ind = 1
with open(filt_bcs) as f:
	for line in f.readlines():
		bc = line.strip()
		bcumis[''.join(['CB:Z:',bc])] = mdict.create("i32:i32") # mdict.create('str:str', 70, 500) #"i64:i32") #dict() 
out_file = sys.argv[4]



#d_let = {'A':'00','C':'01','G':'10','T':'11'}
d_let = {'A':'0','C':'1','G':'2','T':'3'}

#bcumis = defaultdict(lambda: defaultdict(lambda: dict() ))
#bcumis = dict()
def add2bcumis_not(bc, umi, value):
	global bcumis
	#umi1=umi[5]

	#umi1=sys.intern(umi[5])

	#umi1=sys.intern(umi[5:7])

	#umi1=sys.intern(umi[5:9])
	#umi2=sys.intern(umi[10:14])
	#umi3=sys.intern(umi[15:19])

	#umi1=int(''.join([d_let[nt] for nt in list(umi[5:9])]),2)
	#umi2=int(''.join([d_let[nt] for nt in list(umi[10:14])]),2)
	#umi3=int(''.join([d_let[nt] for nt in list(umi[15:19])]),2)

	#umicode=(umi1,umi2,umi3)

	#umicode=umi1*1000000+umi2*1000+umi3

	umicode=int(''.join([d_let[nt] for nt in umi[5::]]))

	#umi1=int(umicode[0:2])
	#umicode=int(umicode)

	if umi1 not in bcumis[bc]:
		bcumis[bc][umi1] = {umicode:value}
	else:
		bcumis[bc][umi1][umicode] = value

	#if umicode not in bcumis[bc]:
	#	bcumis[bc][umicode] = value
	#else:
	#	bcumis[bc][umicode] = False

	#if bc not in bcumis:
	#	bcumis[bc] = {umi1:{umi2:{umi3:value}}}
	#if umi1 not in bcumis[bc]:
	#	bcumis[bc][umi1] = {umi2:{umi3:value}}
	#elif umi2 not in bcumis[bc][umi1]:
	#	bcumis[bc][umi1][umi2] = {umi3:value}
	#else:
	#	bcumis[bc][umi1][umi2][umi3] = False
	
	#elif umi1 not in bcumis[bc]:
	#	bcumis[bc][umi1] = {umi:value}
	#else:
	#	bcumis[bc][umi1][umi] = value


def bad2bcumis(bc, umi):
	global bcumis
	bcumis[bc][umi] = 0 # 0 False ''

def good2bcumis(bc, umi, value):
	global bcumis, keeper, keeper_ind
	bcumis[bc][umi] = keeper_ind # keeper_ind value
	#keeper.append(value)
	keeper.write(value)
	keeper_ind += 1

out = dict()
#out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [int(), str()] ))) #list() # [int(), str()] #, tuple() 
#out_def=defaultdict(lambda: defaultdict(lambda: [int(), str()] ))
#out = defaultdict(lambda: out_def)

chrom_old=''
#d_let = {'A':'0','C':'1','G':'2','T':'3'}
l='ACGT'
umi2code_d = dict(zip([''.join([a,b,c,d])
    for a in l for b in l for c in l for d in l],
    range(256)))
#from textwrap import wrap
#from functools import reduce
#a1 = [l[0:4], l[4:8],l[8:12] ] # wrap(l,4)
#a2 = map(lambda x: umi2code_d[x], a1)
#reduce(lambda x, y: x * 1000 + y, a2)


for line in sys.stdin:
	line = line.strip().split() #"\t" #[5::]
	#strand, chrom, start,cigar,bc,umi = line.strip().split(' ')
	#if (line[4]!="255" or int(line[1])>=256 or line[-2] not in filt_bcs_d or line[-2]=="CB:Z:-" or line[-1]=="UB:Z:-" or line[2]=="chrM"):
	if (line[-2] not in bcumis or line[2]=="chrM"):
		#if (bc not in bcumis or chrom=="chrM"):
		continue

	bc = line[-2] #[5::]
	umi= line[-1] #int(umi)
	#umi = hash(line[-1]) # #[5::] # export PYTHONHASHSEED=1 in terminal before running script
	# way faster but i dont wanna collisions :( 
	# even tho the probability of it is slow, but 
	# i don't wanna risk it
	#umi = int(''.join([d_let[nt] for nt in umi[5::]])) # #[5::] # export PYTHONHASHSEED=1 in terminal before running script
	# found faster
	umi = umi2code_d[umi[5:9]]*1000000 + umi2code_d[umi[9:13]]*1000 + umi2code_d[umi[13:17]] 
	#continue

	#add2bcumis(bc, umi, False)
	#bad2bcumis(bc, umi)
	#continue

	if umi in bcumis[bc]:
		bad2bcumis(bc, umi)
		continue

	strand = line[1][0] 
	if (bc in out) and (umi in out[bc]) and (strand not in out[bc][umi]): 
		bad2bcumis(bc, umi) 
		continue

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

#exit()

if False:
	total_dict=0
	total_umistrings=0

	total_dict += sys.getsizeof(bcumis)
	#print("bcumis",sys.getsizeof(bcumis))
	for bc in bcumis:
		#print(bc, sys.getsizeof(bcumis[bc]))
		#total_dict=0
		#total_umistrings=0
		total_umistrings += sys.getsizeof(bc)
		total_dict += sys.getsizeof(bcumis[bc])
		for umi1 in bcumis[bc]:
			#print(bc, umi1, sys.getsizeof(bcumis[bc][umi1]))
			total_umistrings += sys.getsizeof(umi1)
			#total_dict += sys.getsizeof(bcumis[bc][umi1])
			#for umicode in bcumis[bc][umi1]:
			#	total_umistrings += sys.getsizeof(umicode)
	print("total_dict", total_dict )
	print("total_umistrings", total_umistrings )

	# import psutil
	# pid = os.getpid()
	# memUse = psutil.Process(pid).memory_info()[0]/2.**10 #in kbytes **20 MB **30 GB

	#exit()

proc(chrom)

keeper.seek(0)
result = keeper.read().split(';')
with open(out_file, 'w') as f:
	for bc in bcumis:
		for umi in bcumis[bc]:
			if bcumis[bc][umi]:
				#print(bc, umi, bcumis[bc][umi], file=f)
				print(bc, umi, result[bcumis[bc][umi]], file=f)
exit()