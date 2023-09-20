## covv output from chunks to [starts_out, ends_out]
## cov for BC => from creating additional heavy list to reorganizing what we have
## output bc umi chr? strand? positions which went through filters, or better chunks

#### SAVE COUNTED CIGAR AND TRANSMIT TO THE NEXT STEP!!!
import sys
from collections import defaultdict
import re

def cigar2pos(start_cigar): 
	start, cigar = start_cigar
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
	#starts_out = list()
	#ends_out = list()
	while startss:
		if starts < ends:
			cov += 1
			if cov >= th and not chunk: #switch_cov5:
				chunk = [starts] # True #(starts,)
				#starts_out.append(starts)
			starts = startss.pop()
		elif starts > ends:
			cov -= 1
			if chunk and cov < th : #switch_cov5
				chunk.append(ends) #+= (ends,)
				#ends_out.append(ends)
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
				chunk = [starts] # True #(starts,)
			starts=endss[0]
		cov -= 1
		if chunk and cov < th : 
			chunk.append(ends) # += (ends,)
			#ends_out.append(ends)
			cov5.append(chunk) # += (chunk, )
			chunk = False
		ends = endss.pop()
	#if starts_out:
	#	return([starts_out, ends_out]) #cov5)
	#else:
	#	return(False)
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

chroms = dict()
chrom_ind = 1

bcumis = defaultdict(lambda: dict() )
out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: list()))) #, tuple()
chrom_old=''

for line in sys.stdin:
	line = line.strip().split() #"\t"
	if (line[-2]!="CB:Z:-" and line[-1]!="UB:Z:-" and line[4]=="255" and line[2]!="chrM" and int(line[1])<256):
		bc = line[-2][5::] 
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
				out[bc][umi][strand].append( (start, cigar) ) # += 1 
			elif chrom_old:
				for _bc in out:
					for _umi in out[_bc]:
						if len(out[_bc][_umi]) > 1:
							continue
							# if there are 5+ 5- and 5 of next chr, will pass
							# otherwise bcumis[_bc][_umi] = False needed
						_strand = list(out[_bc][_umi].keys())[0]
						#for _strand in out[_bc][_umi]:
						if len(out[_bc][_umi][_strand]) >= dup_min:
							#print(_bc, _umi)
							# stuff below works for strand check as well
							# as long as there are only 0 and 16 flags
							# 0 forward, 16 reverse
							# could there be other flags with count>th?
							#if (_bc, _umi) not in bcumis:
							if (_bc not in bcumis) or (_umi not in bcumis[_bc]):
								startss = list() 
								endss = list()
								for start_cigar in out[_bc][_umi][_strand]:
									starts, ends = cigar2pos(start_cigar)
									# in c++ better save not absolute but relative nrs, 
									# since size of uint8 << uint32
									# BUT HERE WE ALSO COULD INTRODUCE IT VIA NUMPY
									# not much help, np.uint8 25bytes, np.uint32 28bytes
									startss.extend(starts) #+= (*starts,)
									endss.extend(ends)
								chunks = covv(startss, endss, th = dup_min-1 )
								if chunks:
									bcumis[_bc][_umi] = (_strand, chunks)
							else:
								#bcumis[_bc][_umi] = False
								#thus we could make loops easier but if there are 5+ 5- and 5 of next chr, will pass
								del bcumis[_bc][_umi]

					if len(bcumis[_bc])<cov_min:
						del bcumis[_bc]
					else:
						startss = defaultdict(lambda: list())
						endss = defaultdict(lambda: list())
						for _umi in  bcumis[_bc]:
							_strand, chunks = bcumis[_bc][_umi] #starts_ends
							for chunk in chunks:
								startss[_strand].append(chunk[0])
								endss[_strand].append(chunk[1])
							#starts, ends = starts_ends
							#startss[_strand].extend(starts)
							#endss[_strand].extend(ends)

						starts_ends = dict() #chunks
						for _strand in startss:
							chunkss = covv(startss[_strand], endss[_strand], th=cov_min)
							starts_ends[_strand] = zip(*chunkss)

					for _umi in bcumis[_bc]: #out[_bc]
						_strand, umi_chunks = bcumis[_bc][_umi] #list(out[_bc][_umi].keys())[0]
						uc_starts, uc_ends = zip(*umi_chunks)
						uc_starts.extend(starts_ends[_strand][0])
						uc_ends.extend(starts_ends[_strand][1])
						coords = covv(uc_starts, uc_ends, th=1)
						if coords:
							bcumis[_bc][_umi] = (chrom_old, _strand, coords)
						else:
							bcumis[_bc][_umi] = False

				out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: list())))
				out[bc][umi][strand].append( (start, cigar) ) #  += 1
				chrom_old = chrom #; print('!=,+')
			else: # for the first round (chrom_old == False)
				out[bc][umi][strand].append( (start, cigar) ) #  += 1
				chrom_old = chrom #; print('!=,-')

for _bc in list(out.keys()):
	for _umi in out[_bc]:
		for _strand in out[_bc][_umi]:
			if len(out[_bc][_umi][_strand]) >= dup_min:
				if (_bc, _umi) not in bcumis:
					bcumis[(_bc, _umi)] = (chrom_old, _strand)
				else:
					bcumis[(_bc, _umi)] = False


for bc in bcumis:
	for umi in bcumis[bc]:
		if bcumis[bc][umi];
			chrom, strand, coords = bcumis[bc][umi]
			print(bc, chrom, strand, umi, *coords )