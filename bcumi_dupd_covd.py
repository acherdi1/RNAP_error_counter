## CODE IS TO BE REVISED IF THERE ARE > FLAGS(<256) THAN 0 AND 16 

# TAKES UP ~6GB on AL4_yng while the most in original code is 2.7GB
# 1191s
# maybe separate duplicate filtering again?
# by samtools view | cut -f -2,-1 | uniq -c | >=5
# this way would probably be too ram consuming

# filtering min per bc or per bc_strand??

## output bc umi chr? strand? positions which went through filters, or better chunks
#### SAVE COUNTED CIGAR AND TRANSMIT TO THE NEXT STEP!!!

import sys
from collections import defaultdict
import re
import bisect

# profiling
#from time import time
#https://stackoverflow.com/questions/276052/how-to-get-current-cpu-and-ram-usage-in-python
#import os
#import psutil
#pid = os.getpid()
#memUse = psutil.Process(pid).memory_info()[0]/2.**10 #in kbytes **20 MB **30 GB
#times = defultdict(lambda: list())
#mems = defultdict(lambda: list())



#python -m timeit --setup="x='SRR21398053.324550422 16 chr1 3000096 255 30M238565N120M * 0 0  GCTTTTTTTTTTTTTTTTTTTTTTTTTTGGCTGAGTAGTACTCCATTGTGAAGATGTACCACATTTTCTGTATCCATTCCTCTGTTGAGGGGCATCTGGGTTCTTTCCAGCTTCTGGCTATTATAAATAAGGCTGCTATGAACATAGTGG ,,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF GX:Z:ENSMUSG00000051951 GN:Z:Xkr4 CB:Z:CTGCCATAGGGTTGCA UB:Z:-'" "x.split()[-3]"
#500000 loops, best of 5: 488 nsec per loop
#python -m timeit --setup="x='SRR21398053.324550422 16 chr1 3000096 255 30M238565N120M * 0 0  GCTTTTTTTTTTTTTTTTTTTTTTTTTTGGCTGAGTAGTACTCCATTGTGAAGATGTACCACATTTTCTGTATCCATTCCTCTGTTGAGGGGCATCTGGGTTCTTTCCAGCTTCTGGCTATTATAAATAAGGCTGCTATGAACATAGTGG ,,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF GX:Z:ENSMUSG00000051951 GN:Z:Xkr4 CB:Z:CTGCCATAGGGTTGCA UB:Z:-'" "y=x.rpartition('GN:Z:')" "if y[1]: g=y[2].partition(' ')[0]"
#1000000 loops, best of 5: 215 nsec per loop

#head -100k al4_yng.bam | grep "93M" | wc -l
#78900
#head -100k al4_yng.bam | grep -v "GN:Z:" | wc -l
#33824
#head -100k al4_yng.bam | grep -v 255 | wc -l
#11608
#head -100k al4_yng.bam | grep "CB:Z:-" | wc -l
#2255
#head -100k al4_yng.bam | grep "UB:Z:-" | wc -l
#2733

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
#head -100k al4_yng.bam | grep "93M" | wc -l
#78900
#timeit cigar2pos.py 1_93M 0.3s
#timeit cigar2pos.py 1_63M30S 1.1s
# after introducing if cigar=="93M": return([pos], [pos+93])
# total time decreased only for 1s, meaning that
#this function is not the limiting step at all

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

def proc(chrom):
	global out #bcumis,
	#if 'ACTGTCCAGTTACGAA' in out:
	#	if 'TCTTATTAAATA' in out['ACTGTCCAGTTACGAA']:
	#		print('ACTGTCCAGTTACGAA TCTTATTAAATA', chrom)
	#		print(out['ACTGTCCAGTTACGAA']['TCTTATTAAATA'])
	#if 'ATCACGAGTCATTGCA' in out:
	#	if 'TTATATCCATAC' in out['ATCACGAGTCATTGCA']:
	#		print('ATCACGAGTCATTGCA TTATATCCATAC', chrom)
	#		print(out['ATCACGAGTCATTGCA']['TTATATCCATAC'])
	for bc in out: #range(bc_ind): #
		#if not out[bc]:
		#	continue
		for umi in list(out[bc].keys()): # changed size bc of +out[bc][chrom]
			#if isinstance(out[bc][umi], bool) or isinstance(out[bc][umi], str):
			#	continue
			#if bc in skipped and umi in skipped[bc]:
			#	print('starting', bc, umi, chrom)
			if len(out[bc][umi]) > 1:
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, chrom, "L127 multistrand")
				#	print(out[bc][umi])
				bcumis[bc][umi] = False 
				#out[bc][umi] = False 
				continue
			strand = list(out[bc][umi].keys())[0]
			#if len(out[bc][umi][strand]) < dup_min:
			if out[bc][umi][strand][0] < dup_min:
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, chrom, "L136 not enough dups")
				#	print(out[bc][umi])
				bcumis[bc][umi] = False
				#out[bc][umi] = False
				continue
			#if (bc not in bcumis) or (umi not in bcumis[bc]):
			starts_ends_dups = list()
			start_cigars = out[bc][umi][strand][1]
			start_cigars = start_cigars.split(',')#[1::] ONLY IF out defaultdict str(), ','.joni(out, first start_cigar) => ",.." => empty first, but we switched to normal dict!
			#if bc in skipped and umi in skipped[bc]:
			#	print(bc, umi, 'out[bc][umi][strand][1]', out[bc][umi][strand][1])
			#	print(bc, umi, 'start_cigars', start_cigars)
			for start_cigar in start_cigars:
				starts_ends_dup = cigar2pos(start_cigar)
				#if bc in skipped and umi in skipped[bc]:
				#	print('start_cigar', start_cigar)
				#	print('starts_ends_dup', starts_ends_dup)
				starts_ends_dups.append(starts_ends_dup)
			#if bc in skipped and umi in skipped[bc]:
			#	print('starts_ends_dups', starts_ends_dups)
			#print(starts_ends_dups)
			#exit()
			#print("dups2umi")
			starts_ends_umi = íntersect(starts_ends_dups, th = dup_min ) #-1
			#if bc in skipped and umi in skipped[bc]:
			#	print('starts_ends_umi', starts_ends_umi)
			#print('out', starts_ends_umi)
			if starts_ends_umi:
				#bcumis[bc][chrom][umi] = (strand, starts_ends_umi)
				if chrom in out[bc]:
					out[bc][chrom][umi] = (strand, starts_ends_umi)
				else:
					out[bc][chrom] = {umi:(strand, starts_ends_umi)}
				#if bc in skipped and umi in skipped[bc]:
				#	print('out[bc][chrom][umi]', out[bc][chrom][umi])
			else:
				#out[bc][umi] = False
				#print("122 false") ; print(umi, starts_ends_dups )
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, chrom, "L164 not enough dups in intersect")
				#	#print(out[bc][umi])
				#	#print(starts_ends_dups)
				#	#print(starts_ends_umi)
				bcumis[bc][umi] = False

		if chrom not in out[bc]:
			continue
		#if len(bcumis[bc][chrom])<cov_min:
		if len(out[bc][chrom]) < cov_min :
			#for umi in  bcumis[bc][chrom]:
			#if bc in skipped and chrom=="chrX":
			#	print('len(out[bc][chrom])', len(out[bc][chrom]))
			#	print('out[bc][chrom]', out[bc][chrom])
			#	print(bc, "L176 BC does not have enough UMIs")
			for umi in  out[bc][chrom]:
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, chrom, "L176 BC does not have enough UMIs")
				#	#print(out[bc][chrom])
				bcumis[bc][umi] = False
				#print("131 false")
				#out[bc][umi] = False
			del out[bc][chrom]
			continue

		#len(bcumis[bc][umi])==3
		# that's horrific, 
		#after each chr it would check ALL the umis of each BC present,
		# instead of only working with the output from last chr
		# FIXED !! intermediate bcumis[bc][chrom]

		starts_ends_umis = defaultdict(lambda: list())
		#for umi in  bcumis[bc][chrom]:
		for umi in  out[bc][chrom]:
			#if not bcumis[bc][umi] or len(bcumis[bc][umi])==3:
			#	continue
			#strand, starts_ends_umi = bcumis[bc][chrom][umi]
			strand, starts_ends_umi = out[bc][chrom][umi]
			starts_ends_umis[strand].append(starts_ends_umi)
		#if bc in skipped and chrom=="chrX":
		#	print('starts_ends_umis', starts_ends_umis)
		starts_ends_chrbc = dict()
		for strand in starts_ends_umis:
			starts_ends_chrbc[strand] = íntersect(starts_ends_umis[strand], th=cov_min) #-1
			#if bc in skipped and chrom=="chrX":
			#	print('starts_ends_chrbc[strand]', strand, starts_ends_chrbc[strand])
		#for umi in list(bcumis[bc][chrom].keys()): 
		for umi in list(out[bc][chrom].keys()): 
			#if not bcumis[bc][umi]  or len(bcumis[bc][umi])==3:
			#	continue
			#strand, starts_ends_umi = bcumis[bc][chrom][umi]
			strand, starts_ends_umi = out[bc][chrom][umi]
			#if bc in skipped and umi in skipped[bc]:
			#	print('umi, strand, starts_ends_umi', umi, strand, starts_ends_umi)
			
			if not starts_ends_chrbc[strand]:
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, "L211 BC does not have even a single NT covered with enough UMIs")
				#	#print(out[bc][umi])
				bcumis[bc][umi] = False
				#print("163 false") ; print(umi, strand, starts_ends_umis[strand])
				#out[bc][umi] = False
				#del bcumis[bc][chrom][umi]
				del out[bc][chrom][umi]
				continue
			coords = íntersect([starts_ends_chrbc[strand], starts_ends_umi], th=2)
			#if bc in skipped and umi in skipped[bc]:
			#	print('umi, coords', umi, coords)
			if coords:
				#bcumis[bc][chrom][umi] = (strand, coords)
				out[bc][chrom][umi] = (strand, coords)
			else:
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, chrom, "L225 UMI does not participate in any covered enough regions in BC (check borders?)")
				#	#print(starts_ends_umi)
				#	#print(starts_ends_umis[strand])
				#	#print(starts_ends_chrbc[strand])
				bcumis[bc][umi] = False
				#print("174 false") ; print(umi,starts_ends_chrbc[strand], starts_ends_umi )
				#del bcumis[bc][chrom][umi]
				#out[bc][umi] = False
				del out[bc][chrom][umi]


		#if len(bcumis[bc][chrom])<cov_min:
		if len(out[bc][chrom])<cov_min:
			#for umi in  bcumis[bc][chrom]:
			#	bcumis[bc][umi] = False
			#if bc in skipped and chrom=="chrX":
			#	print(bc, umi, chrom, "L240 no UMI enough in BC (could be the case with TH-1 UMIs left)")
			for umi in  out[bc][chrom]:
				#out[bc][umi] = False
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, chrom, "L240 no UMI enough in BC (could be the case with TH-1 UMIs left)")
				bcumis[bc][umi] = False
				#print("187 false")
			del out[bc][chrom]
			continue

		#print(2)
		#for umi in  bcumis[bc][chrom]:
		for umi in  out[bc][chrom]:
			#strand, coords = bcumis[bc][chrom][umi]
			strand, coords = out[bc][chrom][umi]
			#print(bc, chrom, strand, umi, *coords )
			#bc_ = l_bcs[bc] #; _chr = list(chroms.keys())[chrom]
			#umi_ = ''.join([d_let_inv[i] for i in list(str(umi))])
			#print(bc, chrom, strand, umi, *coords )
			#bcumis[bc][umi] = True #(chrom, strand, coords)
			#out[bc][umi] = True #(chrom, strand, coords)
			#out[bc][umi] = ' '.join([chrom, strand, str(coords[0]), str(coords[1])])
			bcumis[bc][umi] = ' '.join([chrom, strand, str(coords[0]), str(coords[1])])
			#if bc in skipped and umi in skipped[bc]:
			#	print('bcumis[bc][umi]',bcumis[bc][umi])

		#del bcumis[bc][chrom]
		del out[bc][chrom]

def add2out(bc, umi, strand, start_cigar):
	global out
	if bc in out:
		if umi in out[bc]:
			if strand in out[bc][umi]:
				out[bc][umi][strand][0] += 1
				out[bc][umi][strand][1] = ','.join([out[bc][umi][strand][1], start_cigar])
			else:
				#if bc in skipped and umi in skipped[bc]:
				#	print(bc, umi, chrom, "L271 there is another strand already")
				bcumis[bc][umi] = False
				#print("add2out false")
		else:
			out[bc][umi] = {strand:[1, start_cigar]}
	else:
		out[bc] = {umi:{strand:[1, start_cigar]}}


#skipped=dict() #"bcumi_absent"
#with open(sys.argv[5]) as f:
#	for line in f:
#		bc,umi = line.strip().split()
#		if bc not in skipped:
#			skipped[bc] = {umi: 1}
#		else:
#			skipped[bc][umi] = 1


pattern = re.compile(r'(\d+)([MIDNS])') 

#th = int(sys.argv[1])
dup_min = int(sys.argv[1]) #5
cov_min = int(sys.argv[2]) #5

#if False: # if reading from reads.txt
filt_bcs=sys.argv[3]
filt_bcs_d=dict()
with open(filt_bcs) as f:
	bc_ind = 1 #0
	for line in f.readlines():
		bc = line.strip()
		filt_bcs_d[bc] = bc_ind
		#bc_ind += 1
#l_bcs = list(filt_bcs_d.keys())

out_file = sys.argv[4]
#bcumis_excl=sys.argv[4]

#chroms = dict()
#chrom_ind = 0

d_let = {'A':'1','C':'2','G':'3','T':'4'}
d_let_inv = {v:k for k, v in d_let.items()}

bcumis = defaultdict(lambda: defaultdict(lambda: dict() ))
#bcumis = dict()
out = dict()
#out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [int(), str()] ))) #list() # [int(), str()] #, tuple() 
#out = [None]*bc_ind


chrom_old=''


for line in sys.stdin:
	line = line.strip().split() #"\t"
	#if line[2]!="chr1" or line[-2]!="CB:Z:ACACAGTAGAGAGGGC" : #if line[2]!="chr1" or 
	#	continue
	# if reading from reads.txt
	#bc, chr_str, umi, start, cigar, _, _ = line
	#chrom, strand = chr_str.split('_')
	#start_cigar = '_'.join([start, cigar])
	#if False: # if reading from reads.txt
	if (line[4]!="255" or int(line[1])>=256 or line[-2]=="CB:Z:-" or line[-1]=="UB:Z:-" or line[2]=="chrM"):
		continue
	bc = line[-2][5::] 
	if bc not in filt_bcs_d:
		continue

	#bc = filt_bcs_d[bc] # replacing bc with its index (string to int)
	#if not out[bc]:
	#	out[bc] = defaultdict(lambda: defaultdict(lambda: list()))
	umi = line[-1][5::] #int(''.join([d_let[i] for i in list(line[-1][5::])]))
	#if ((bc in out) and (umi in out[bc]) and
	#	(isinstance(out[bc][umi], bool) or isinstance(out[bc][umi], str))):
	if (bc in bcumis) and (umi in bcumis[bc]):
		#out[bc][umi] = False
		#if bc in skipped and umi in skipped[bc]:
		#	print(bc, umi, chrom, "L349 already was in another chrom (here could be wrong chrom annotated, from prev loop)")
		bcumis[bc][umi] = False
		continue
	strand = line[1][0] # 0 forward, 16 reverse; could there be other flags with count>th?
	
	chrom = line[2] # below replacing chrom with its index
	#if chrom not in chroms:
	#	chroms[chrom] = chrom_ind
	#	chrom_ind += 1
	#chrom = chroms[chrom]

	#start = int(line[3]); cigar = line[5] #int(
	start_cigar = '_'.join([line[3], line[5]])

	if chrom == chrom_old:
		#out[bc][umi][strand].append( start_cigar) # (start, cigar) ) #
		add2out(bc, umi, strand, start_cigar)

	elif chrom_old:
		proc(chrom_old)
		#if chrom == "chr3":
		#	exit()
		#out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [int(), str()] ))) #list() [int(), str()]
		#out = [None]*bc_ind
		#out[bc] = defaultdict(lambda: defaultdict(lambda: list()))
		#out[bc][umi][strand].append( start_cigar) # (start, cigar) ) #
		
		#out[bc][umi][strand][0] += 1
		#out[bc][umi][strand][1] = ','.join([out[bc][umi][strand][1], start_cigar])
		out = {bc:{umi:{strand:[1, start_cigar]}}}
		#add2out(bc, umi, strand, start_cigar)
		chrom_old = chrom 

	else: # for the first round (chrom_old == False)
		#out[bc][umi][strand].append( start_cigar) # (start, cigar) ) #
		#out[bc][umi][strand][0] += 1
		#out[bc][umi][strand][1] = ','.join([out[bc][umi][strand][1], start_cigar])
		add2out(bc, umi, strand, start_cigar)
		chrom_old = chrom 

proc(chrom)

#with open(bcumis_excl, 'w') as f:
#	#for bc in bcumis:
#	#	for umi in bcumis[bc]:
#	#		if not bcumis[bc][umi]:
#	for bc in out:
#		for umi in out[bc]:
#			if not out[bc][umi]:
#				#_bc = l_bcs[bc]
#				#_umi = ''.join([d_let_inv[i] for i in list(str(umi))])
#				print(bc, umi, file=f)


#for bc in out:
#	for umi in out[bc]:
#		if out[bc][umi]:
#			print(bc, umi, out[bc][umi])

#exit()
with open(out_file, 'w') as f:
	for bc in bcumis:
		for umi in bcumis[bc]:
			if bcumis[bc][umi]:
				print(bc, umi, bcumis[bc][umi], file=f)
exit()

#for bc in out:
#	for umi in bcumis[bc]:
#		if bcumis[bc][umi]:
#			chrom, strand, coords = bcumis[bc][umi]
#			#_bc = l_bcs[bc] #; _chr = l_chrs[chrom]
#			#_umi = ''.join([d_let_inv[i] for i in list(str(umi))])
#			print(bc, chrom, strand, umi, *coords )