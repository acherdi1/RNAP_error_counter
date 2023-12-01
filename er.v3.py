# memory improvement: save BC as indices as in bcumi_dupd_cov.py


# when len_cov is dict and counted out of sript:
# time 458s
# RAM 2.43GB

# cov to str => RAM 2.42 , 
#but comparison in strings one by one letter?
# time 462s, neglectable increase


#from collections import defaultdict
import sys
import re
#from microdict import mdict
#import numpy as np

max_chr_len = 195154279
dup_min = int(sys.argv[1]) #4
cov_min = int(sys.argv[2]) #5
vars_exp = int(sys.argv[3]) #2
#work_dir = sys.argv[1]
#input_file = sys.argv[4] #'/'.join([work_dir, "reads.txt"])
#bcumis_file= sys.argv[5] #'/'.join([work_dir,"bcumi_covd.txt"])
cov_file = sys.argv[4] #'interm4'
address_subs= sys.argv[5] #'/'.join([work_dir, "subs.txt"])
address_er= sys.argv[6] #'/'.join([work_dir, "er.txt"])

with open(address_subs, 'w') as f:
	pass

out = dict()

#chroms = dict()
#chrom_ind = 0

# with open(cov_file) as f:
# 	for line in f: 
# 		line = line.strip()
# 		info, starts, ends = line.split("[")
# 		#bc, chrom, strand, umi = info.split()
# 		bc, umi, chrom, strand = info.split()
# 		strand = int(strand[0])
# 		umi = int(umi) + max_chr_len
# 		#if chrom not in chroms:
# 		#	chroms[chrom] = chrom_ind
# 		#	chrom_ind += 1
# 		#chrom = chroms[chrom]
# 		starts = [int(i) for i in starts.strip('] ').split(', ')]
# 		ends = [int(i) for i in ends.strip('] ').split(', ')]
# 		if chrom not in out:
# 			out[chrom] = dict()
# 		if (bc,strand) not in out[chrom]:
# 			out[chrom][bc,strand] = {umi:(starts, ends)}
# 		else:
# 			out[chrom][bc,strand][umi] = (starts, ends)



# def unpack_startsends(chrom,bc_strand,umi):
# 	global out
# 	starts, ends = out[chrom][bc_strand][umi]
# 	out[chrom][bc_strand][umi] = dict()
# 	#del out[chrom][bc_strand][umi]
# 	#out[chrom][bc_strand][umi] = mdict.create("str:str", key_len=None, val_len=None)
# 	for chunk in range(len(starts)):
# 		for pos in range(starts[chunk], ends[chunk]):
# 			#out[chrom][bc_strand][umi][pos] = [0,0,0,0,0] #{1: 0, 2: 0, 3: 0, 4: 0, 0: 0} #[nt[i]] = 1
# 			out[chrom][bc_strand][umi][pos] = b""

def unpack_startsends(chrom): #,bc_strand,umi):
	global out
	with open(cov_file) as f:
		chrom_check = chrom + " "
		for line in f: 
			if chrom_check not in line:
				continue
			line = line.strip()
			info, starts, ends = line.split("[")
			bc, umi, chrom, strand = info.split()
			strand = int(strand[0])
			umi = int(umi) + max_chr_len
			starts = [int(i) for i in starts.strip('] ').split(', ')]
			ends = [int(i) for i in ends.strip('] ').split(', ')]
			if chrom not in out:
				out[chrom] = {(bc,strand):{umi:[(starts, ends), b""]}} # dict()
			elif (bc,strand) not in out[chrom]:
				out[chrom][bc,strand] = {umi:[(starts, ends), b""]} #(starts, ends)
			else:
				out[chrom][bc,strand][umi] = [(starts, ends), b""] #(starts, ends)

			#for chunk in range(len(starts)):
			#	for pos in range(starts[chunk], ends[chunk]):
			#		out[chrom][bc,strand][umi][pos] = b""




def pileup(pos, cigar, seq):
	global out
	pos_in = 0
	cigar_ops = pattern.findall(cigar) 
	if cigar_ops[0][1] == 'S':
		pos_in += int(cigar_ops[0][0])
		cigar_ops = cigar_ops[1:]
	for length, op in cigar_ops:
		length = int(length)
		if op == 'M':
			#coords[pos_in:(pos_in+length)] = range(pos, pos + length)
			#out[chrom][bc,strand][umi][1] += bytes("".join([pos, "-", pos+length, ","]), encoding='utf-8')
			#out[chrom][bc,strand][umi][2] += bytes("".join([seq[pos_in:pos_in+i], ","]), encoding='utf-8')
			coord = [pos, pos+length]
			starts, ends = out[chrom][bc,strand][umi][0]

			for se_ind in range(len(starts)):
				#print(se_ind)
				#print(starts[se_ind], ends[se_ind], coord, seq)
				if coord[1] < starts[se_ind]:
					continue
				if coord[0] > ends[se_ind]:
					continue #break

				#print(starts[se_ind], ends[se_ind], coord, seq)
				start_in = starts[se_ind] - coord[0]
				if start_in < 0:
					start = coord[0]
					start_in = pos_in #0
				else:
					start = starts[se_ind]
					start_in += pos_in

				end_in = ends[se_ind] - coord[1]
				if end_in > 0:
					end_in = coord[1]-coord[0]+1 +pos_in
					end = coord[1] +1
				else:
					end_in += coord[1]-coord[0]+1 +pos_in
					end = ends[se_ind] +1

				#length_out = end - start
				seq_out = seq[start_in:end_in]
				#print(start,length_out, seq_out)
				out[chrom][bc,strand][umi][1] += bytes("".join([str(start), ":", seq_out, ","]), encoding='utf-8') # "-", str(length), 
				#pos_in += length

			#out[chrom][bc,strand][umi][1] += bytes("".join([str(pos), ":", seq[pos_in:pos_in+length], ","]), encoding='utf-8') # "-", str(length), 
			
			#for i in range(length):
			#	if pos+i in out[chrom][bc,strand][umi]:
			#		#nt_ind = nt2int[seq[pos_in+i]]
			#		#out[chrom][bc,strand][umi][pos+i][nt_ind] += 1
			#		out[chrom][bc,strand][umi][pos+i] += bytes(seq[pos_in+i], encoding='utf-8') #= "".join([out[chrom][bc,strand][umi][pos+i], seq[pos_in+i]])
			
			pos_in += length
		elif op == 'I':
			pos_in += length
			continue
		pos += length
	out[chrom][bc,strand][umi][1] += b" "
	#out[chrom][bc,strand][umi][2] += b" "

def leave_only_max(d):
	#print(d)
	dmax=[0,0]
	for i in range(5):
		if d[i]>dmax[1]:
			dmax=[i,d[i]]
	#print(dmax)
	return(dmax)

def umi_consensus(bc_strand, chrom, umi): #, dup_min
	global out
	for pos in list(out[chrom][bc_strand][umi].keys()):
		#print(chrom, bc_strand, umi, pos)
		nt, cov = leave_only_max(out[chrom][bc_strand][umi][pos])
		#if cov < dup_min or not nt:
		#	del out[chrom][bc_strand][umi][pos]
		#else:
		#	out[chrom][bc_strand][umi][pos] = nt
		if cov >= dup_min and nt:
			if pos in out[chrom][bc_strand]:
				out[chrom][bc_strand][pos][umi] = nt
			else:
				out[chrom][bc_strand][pos] = {umi: nt}
		del out[chrom][bc_strand][umi][pos]

	#if not out[chrom][bc_strand][umi]:
	del out[chrom][bc_strand][umi]

def variance(bc_strand, chrom): 
	global out, len_cov
	len_cov_chr = 0
	len_cov_err_chr = 0
	for pos in list(out[chrom][bc_strand].keys()):
		if len(out[chrom][bc_strand][pos]) < cov_min:
			del out[chrom][bc_strand][pos]
			continue
		len_cov_chr += 1
		d = dict() 
		for umi in list(out[chrom][bc_strand][pos].keys()):
			if out[chrom][bc_strand][pos][umi] in d:
				d[out[chrom][bc_strand][pos][umi]] += 1
			else:
				d[out[chrom][bc_strand][pos][umi]] = 1
		variants = len(d)
		if variants != vars_exp or min(d.values())>1:
			del out[chrom][bc_strand][pos]
			continue
			# 1) remove those where mismatch is repeated more than once on the same position!
			# so when more than one umi contains mismatch on the same position
			# bc pol2 error hardly could happen twice at the same position
			# meaning that mismatch probably happened due to some other reasons
			# 2) remove if there are more than 2 variants per min position
			# same reasoning
		len_cov_err_chr += 1
		out[chrom][bc_strand][pos] = sorted(d, key=d.get, reverse=True) #d
		out[chrom][bc_strand][pos] = [int2nt[i] for i in out[chrom][bc_strand][pos]]
		# mismatch from which to which nt
		#[1] " AC" " TA"

	#return(len_cov_chr, len_cov_err_chr)
	_bc = bc_strand[0]
	if _bc in len_cov:
		len_cov[_bc][0] += len_cov_err_chr ; len_cov[_bc][1] += len_cov_chr
	else:
		len_cov[_bc] = [len_cov_err_chr, len_cov_chr]

def record_substitutions(address, chrom):
	with open(address, "a") as f:
		for bc_strand in out[chrom]:
			for pos in out[chrom][bc_strand]:
				print(*bc_strand, chrom, pos, *out[chrom][bc_strand][pos], file = f)

def record_errorrate(address): #er, 
	with open(address, "w") as f:
		for bc in len_cov:
			if len_cov[bc][1]:
				print(bc, len_cov[bc][0]/len_cov[bc][1], *len_cov[bc], file = f) #errorRate


def across_chrom(chrom):
	global out
	#for chrom in out:
	print("chrom",chrom)
	for bc_strand in out[chrom]:
		print("bc strand",bc_strand)
		for umi in list(out[chrom][bc_strand].keys()):
			print("umi",umi)
			seqs = str(out[chrom][bc_strand][umi][1], encoding="utf-8")
			seqs = [i.split(":") for i in seqs.split(",")][:-1]
			seqs = [[int(i[0]),i[1]] for i in seqs]
			interm_dict = dict()
			print(interm_dict)
			#print("seqs",seqs)
			for s in seqs:
				if s[0] not in interm_dict:
					interm_dict[s[0]] = [s[1]]
				else:
					interm_dict[int(s[0])].append(s[1])
					#for pos_in in range(len(s[1])):
					#	if pos_in == len(interm_dict[s[0]]): # if current seq longer than prev
					#		interm_dict[s[0]] = "".join([interm_dict[s[0]], s[1][pos_in:]])
					#		break
					#	if s[1][pos_in]==interm_dict[s[0]][pos_in]:
					#		continue
					#	else:
			#interm_dict[103826667] = ["AGTG", "AGAGGT"]
			#interm_dict[103826670] = ["GGTACTTGTGAGC"]
			for i in interm_dict:
				interm_dict[i].append(0)
			#print("interm_dict",interm_dict)
			while True:
				start=min(interm_dict.keys()) # min - for umis, not duplicates, or better do it on the next step?
				pos = start
				#pos_in = 0
				seqs = [interm_dict[pos]]
				interm_dict.pop(pos)
				#print(interm_dict)
				ref = ""
				#print("pos seqs",pos, seqs)
				here = [ i[k[-1]] for k in seqs for i in k[:-1]]
				#print(here)

				if len(set(here)) == 1:
					ref += here[0]
				else:
					b = {}
					for item in here:
						b[item] = b.get(item, 0) + 1
					#print(b)
					ref=''.join([ref,"(",*[str(v)+k for k,v in b.items()], ")" ])
				#print(ref)

				for seq in seqs:
					seq[-1]+=1
				while True:
					pos += 1
					if pos in interm_dict:
						seqs.append(interm_dict[pos])
						interm_dict.pop(pos)
					#print(interm_dict)
					here = [i[k[-1]] for k in seqs for i in k[:-1] if len(i)>k[-1]]
					#print(here, set(here))
					if not here:
						break
					if len(set(here)) == 1:
						ref += here[0]
					else:
						b = {}
						for item in here:
							b[item] = b.get(item, 0) + 1
						#print(b)
						ref=''.join([ref,"(",*[str(v)+k for k,v in b.items()], ")" ])
					#print(ref)
					for seq in seqs:
						seq[-1]+=1
					#print(seqs)
				#print(ref)
				# out[chrom][bc_strand][start] = out[chrom][bc_strand].get(start,[]).append(ref)
				if start not in out[chrom][bc_strand]:
					out[chrom][bc_strand][start] = [ref]
					#if start==103826666:
					#	out[chrom][bc_strand][start].append('TAGTGG(4T1A)ACTTGTG')
					#	out[chrom][bc_strand][start].append('TAGTGG(4T1A)ACTTGTG')
					#	out[chrom][bc_strand][103826670]=['GGTACTTGT'] #GAGCCAT
				else:
					out[chrom][bc_strand][start].append(ref)
				#print(interm_dict)
				if not interm_dict:
					del out[chrom][bc_strand][umi]
					break
		#print(out)

		# first compare if the strings are the same on the length of minlen and befor next start

		#g = iter(sorted(out[chrom][bc_strand].keys()))
		#start= next(g) #min(out[chrom][bc_strand])
		for start in sorted(out[chrom][bc_strand].keys()):
			if start not in out[chrom][bc_strand]:
				continue
			pos = start

			seqs = [out[chrom][bc_strand][pos]]
			seqs[-1].append(0)
			#print("seqs2",seqs)
			out[chrom][bc_strand].pop(pos)
			ref = ""
			#here = [ i[k[-1]] for k in seqs for i in k[:-1]]
			here = list()
			for k in seqs: 
				for i in range(len(k[:-1])): 
					if len(k[i])>k[-1]:
						if k[i][k[-1]]!="(":
							here.append(k[i][k[-1]])
						else:
							nt=k[i][k[-1]:].split(")")[0]+")"
							here.append(nt)
							k[i]=k[i][len(nt)-1:]
				k[-1] += 1
			#print(here)
			s_here = set(here)
			if len(s_here) == 1:
				ref = ''.join([ref, str(len(here)),here[0]])
			else:
				b = {}
				for item in here:
					b[item] = b.get(item, 0) + 1
				ref=''.join([ref,"(",*[str(v)+k for k,v in b.items()], ")" ])
			#print(ref)
			#break
			#for seq in seqs:
			#	seq[-1]+=1
			while True:
				pos += 1
				if pos in out[chrom][bc_strand]:
					seqs.append(out[chrom][bc_strand][pos])
					seqs[-1].append(0)
					out[chrom][bc_strand].pop(pos)
				#print(seqs)
				#print(out[chrom][bc_strand])
				here = list()
				for k in seqs: 
					for i in range(len(k[:-1])): 
						if len(k[i])>k[-1]:
							if k[i][k[-1]]!="(":
								here.append(k[i][k[-1]])
							else:
								nt=k[i][k[-1]:].split(")")[0]+")"
								here.append(nt)
								k[i]=k[i][len(nt)-1:]

					k[-1] += 1
				#print(seqs)
				#here = [i[k[-1]] if i[k[-1]]!="(" else i[1+k[-1]:].split(")")[0] ]
				#print(here)
				#print(seqs)
				#break
				#print(here, s_here)
				#break
				if not here:
					break
				s_here = set(here)
				if len(s_here) == 1:
					ref = ''.join([ref, str(len(here)),here[0]])
				else:
					#if any([i.startswith("(") for i in s_here ]):
					#	print("!!!") 
					#	break
					#	print(here.index("("))
					#	#while True:
					#	#	nr_seqs = 
					#	break
					b = {}
					for item in here:
						b[item] = b.get(item, 0) + 1
					#print(b)
					ref=''.join([ref,"(",*[str(v)+k for k,v in b.items()], ")" ])
				#print(ref)
				#for seq in seqs:
				#	seq[-1]+=1
				#print(seqs)
			#break
			out[chrom][bc_strand][start] = ref
			print(start, ref)
		#print(out)


chrom_old=''
pattern = re.compile(r'(\d+)([MIDNS])') 
nt2int = {"N":0, "A":1,"C":2,"T":3,"G":4}
int2nt = {v: k for k, v in nt2int.items()}

len_cov = dict() 

l='ACGT'
umi2code_d = dict(zip([''.join([a,b,c,d])
    for a in l for b in l for c in l for d in l],
    range(256)))
ind = iter(range(16))
for a in l:
	for b in l:
		umi2code_d[''.join([a,b])] = next(ind)


for line in sys.stdin:
	line = line.strip().split() #"\t"
	#if line[2]!="chr1" or line[-2]!="CB:Z:ACACAGTAGAGAGGGC" :
	#	continue


	# if reading from reads.txt
	#bc, chr_str, umi, start, cigar, seq, _ = line
	#chrom, strand = chr_str.split('_')
	#seq = list(seq)
	#start = int(start)
	#strand = int(strand)
	#start_cigar = '_'.join([start, cigar])
	#if (chrom not in out) or ((bc,strand) not in out[chrom]) or (umi not in out[chrom][bc,strand]):
	#	print(umi)
	#	continue
	#if False: # if reading from reads.txt
	if line[-2]=="CB:Z:-" or line[-1]=="UB:Z:-" or line[4]!="255" or line[2]=="chrM" or int(line[1])>=256:
		continue
	chrom = line[2]
	#if chrom not in chroms:
	#	continue
	#chrom = chroms[chrom]
	bc = line[-2] ; umi= line[-1] ; strand = int(line[1][0])
	umi = umi2code_d[umi[5:9]]*1000000 + umi2code_d[umi[9:13]]*1000 + umi2code_d[umi[13:17]] + max_chr_len # max size of chr in mice, so it would not intersect with pos after umi_consensus in out

	if chrom != chrom_old:
		if chrom_old:
			#print('here', len(out[chrom_old]))
			across_chrom(chrom_old)
			#for bc_strand in out[chrom_old]:
			#	#x = out[chrom_old][bc_strand] # would be the same dict, no
			#	for _umi in list(out[chrom_old][bc_strand].keys()):
			#		umi_consensus(bc_strand, chrom_old, _umi) # _bc _chr _umi
			#	variance(bc_strand, chrom_old)
			#record_substitutions(address_subs, chrom_old)
			del out[chrom_old]
		chrom_old = chrom
		#for bc_strand in out[chrom]:
		#	for _umi in out[chrom][bc_strand]:
		unpack_startsends(chrom) #, bc_strand, _umi)
		#exit()

	#if isinstance(out[chrom][bc,strand][umi], tuple):
	#	unpack_startsends(bc,chrom,strand,umi)

	if (chrom not in out) or ((bc,strand) not in out[chrom]) or (umi not in out[chrom][bc,strand]):
		continue
	#print(umi)

	start = int(line[3]); cigar = line[5]
	seq = line[9] # list(line[9]) # [nt2int[i] for i in list(line[9])] 
	# but is it better? strings from different lists are probably 
	# the same objects, meaning that "A" from all its positions 
	# in all seqs would take total 64 bytes

	pileup(start, cigar, seq)

#print(out)

exit()

for bc_strand in out["chr7"]:
	print(bc_strand)
	for umi in out["chr7"][bc_strand]:
		print(umi)
		print(out["chr7"][bc_strand][umi][0])
		print(str(out["chr7"][bc_strand][umi][1], encoding='utf-8'))

exit()	

if chrom in out:
	for bc_strand in out[chrom]:
		for _umi in list(out[chrom][bc_strand].keys()):
			umi_consensus(bc_strand, chrom, _umi) # _bc _chr _umi
		variance(bc_strand, chrom) 

	record_substitutions(address_subs, chrom)

record_errorrate(address_er)
exit()



# AGGCCACAGGTGCTAG 0 1 TGCGATGTTAAG                                        
# {28277446: [0, 0, 0, 0, 0], 28277447: [0, 0, 0, 0, 0], 
# 28277448: [0, 0, 0, 0, 0], 28277449: [0, 0, 0, 0, 0], 
# 28277450: [0, 0, 0, 0, 0], 28277451: [0, 0, 0, 0, 0], 
# 28277452: [0, 0, 0, 0, 0], 28277453: [0, 0, 0, 0, 0]}
# 28277345        91M2S
# 28277361        93M
# 28277422        63M1D30M
# 28277431        93M
# 28277432        93M
# 28277507        93M