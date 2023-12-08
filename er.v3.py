import sys
import re

dup_min = int(sys.argv[1])
cov_min = int(sys.argv[2])
vars_exp = int(sys.argv[3]) 
cov_file = sys.argv[4] 
address_er= sys.argv[5] 

try:
	print_var=sys.argv[6].split(" ")
	print_var[1] = int(print_var[1])
except:
	print_var=0

def unpack_startsends(chrom): 
	global out
	with open(cov_file) as f:
		chrom_check = chrom + " "
		for line in f: 
			if chrom_check not in line:
				continue
			line = line.strip()
			info, starts, ends = line.split("[")
			bc, umi, chrom, strand = info.split()
			if not bc.startswith("CB"):
				bc = ''.join(["CB:Z:",bc])
			strand = int(strand[0])
			if umi[0].isalpha():
				umi = umi2code_d[umi[0:4]]*1000000 + umi2code_d[umi[4:8]]*1000 + umi2code_d[umi[8:12]] + max_chr_len # max size of chr in mice, so it would not intersect with pos after umi_consensus in out
			else:
				umi = int(umi) + max_chr_len
			starts = [int(i) for i in starts.strip('] ').split(', ')]
			ends = [int(i) for i in ends.strip('] ').split(', ')]
			#for i in range(len(starts)):
			#	if starts[i]<=72712169 and ends[i]>=72712169:
			#		print(bc, umi, chrom, strand, starts, ends)
			if chrom not in out:
				out[chrom] = {(bc,strand):{umi:[(starts, ends), b""]}} # dict()
			elif (bc,strand) not in out[chrom]:
				out[chrom][bc,strand] = {umi:[(starts, ends), b""]} #(starts, ends)
			else:
				out[chrom][bc,strand][umi] = [(starts, ends), b""] #(starts, ends)

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
			coord = [pos, pos+length]
			starts, ends = out[chrom][bc,strand][umi][0]

			for se_ind in range(len(starts)):
				if (coord[1] < starts[se_ind]) or (coord[0] > ends[se_ind]):
					continue

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

				seq_out = seq[start_in:end_in]
				out[chrom][bc,strand][umi][1] += bytes("".join([str(start), ":", seq_out, ","]), encoding='utf-8') # "-", str(length), 
			#print("pileup", chrom, bc, umi, length, op, coord, starts, ends, out[chrom][bc,strand][umi][1])
			pos_in += length
		elif op == 'I':
			pos_in += length
			continue
		pos += length
	out[chrom][bc,strand][umi][1] += b" "

def umi2code_fnc(umi):
	umi2 = ''.join([code2umi[(umi-max_chr_len)//1000000], 
		code2umi[(umi-max_chr_len)%1000000//1000],
		code2umi[(umi-max_chr_len)%1000]])
	return umi2

def across_chrom(chrom):
	global out
	for bc_strand in out[chrom]:
		if bc_strand[0] not in cell_cov:
			cell_cov[bc_strand[0]]=[0,0,0]

		for umi in list(out[chrom][bc_strand].keys()):
			dups = str(out[chrom][bc_strand][umi][1], encoding="utf-8").split(",")
			dups = [dup.split(":") for dup in dups][:-1]
			dups = [[int(dup[0]),dup[1]] for dup in dups]

			dups_dict = dict()
			for dup in dups:
				pos, seq = dup
				dups_dict[pos] = dups_dict.get(pos,[]) + [seq] #.append(seq)
			for i in dups_dict:
				dups_dict[i].append(0)
			#print(chrom, bc_strand, umi, "dups", dups, "dups_dict", dups_dict)
			# for each non-interrupted region consisting of multiple duplicates:
			while True: 
				#try:
				#print( 'non-interrupted region', dups_dict )
				start=min(dups_dict.keys()) # start of uninterrupted region
				#except:
				#	break
				ref = ""
				seqs = list()
				pos = start -1

				# while the region is uninterrupted, for each next position
				while True:
					pos += 1
					if pos in dups_dict:
						seqs.append(dups_dict.pop(pos))

					nts_at_pos = [i[k[-1]] for k in seqs for i in k[:-1] if len(i)>k[-1]]
					nt_types_at_pos = set(nts_at_pos)
					if print_var and pos==print_var[1]:
						print(umi2code_fnc(umi), pos, 
							''.join([''.join(s) for s in [(str(nts_at_pos.count(i)),i) for i in nt_types_at_pos]]))				
					if not nts_at_pos:
						break
					if len(nt_types_at_pos) == 1:
						ref += nts_at_pos[0]
					else:
						b = {nt_type:nts_at_pos.count(nt_type) for nt_type in nt_types_at_pos}
						ref=''.join([ref,"(",*[str(v)+k for k,v in b.items()], ")" ])
					for seq in seqs:
						seq[-1]+=1

				out[chrom][bc_strand][start] = out[chrom][bc_strand].get(start,[]) + [ref] #.append(ref)
				if not dups_dict:
					del out[chrom][bc_strand][umi]
					break

		# first compare if the strings are the same on the length of minlen and befor next start
		#g = iter(sorted(out[chrom][bc_strand].keys()))
		#start= next(g) #min(out[chrom][bc_strand])

		for start in sorted(out[chrom][bc_strand].keys()):
			if start not in out[chrom][bc_strand]:
				continue
			pos = start-1
			seqs=[]
			ref = ""

			while True:
				pos += 1
				if pos in out[chrom][bc_strand]:
					seqs.append(out[chrom][bc_strand][pos])
					seqs[-1].append(0) # append counter
					out[chrom][bc_strand].pop(pos)
				nts_at_pos = list()
				for seqs_with_same_start in seqs: 
					pos_in = seqs_with_same_start[-1]
					for seq_nr in range(len(seqs_with_same_start[:-1])): 
						seq = seqs_with_same_start[seq_nr]
						if pos_in >= len(seq):
							continue # here the seq could also be dropped from seqs
						nt = seq[pos_in]
						if seq[pos_in]=="(":
							nt=seq[pos_in:].split(")")[0]+")"
							seqs_with_same_start[seq_nr]=seq[len(nt)-1:]
						nts_at_pos.append(nt)

					seqs_with_same_start[-1] += 1 # pos_in+1 for those seqs with same start pos
				if not nts_at_pos:
					break
				if len(nts_at_pos)<cov_min:
					continue
				nt_types_at_pos = set(nts_at_pos)
				if len(nt_types_at_pos) == 1:
					#ref = ''.join([ref, str(len(nts_at_pos)),nts_at_pos[0]])
					cell_cov[bc_strand[0]][0] += len(nts_at_pos)
					cell_cov[bc_strand[0]][1] += 1
				else:
					b = {nt_type:nts_at_pos.count(nt_type) for nt_type in nt_types_at_pos}
					#ref=''.join([ref,"(",*[str(v)+k for k,v in b.items()], ")" ])
					pcr_ers=list()
					for nt, cnt in list(b.items()):
						if not nt.startswith("("):
							continue
						pcr_vars = pattern_pcr_vars.findall(nt)
						pcr_vars = [(int(i[0]), i[1]) for i in pcr_vars]
						maxx = (0,'')
						for pcr_var in pcr_vars:
							if pcr_var[0]>maxx[0]:
								maxx = pcr_var
						if maxx[0] >= dup_min and ((2*maxx[0]-sum([i[0] for i in pcr_vars]))>=(maxx[0]/2)):
							b[maxx[1]] = b.get(maxx[1],0) + 1 #int(maxx[0])
							pcr_ers.append(''.join([str(cnt),nt]))
						del b[nt]
					#if len(b)==1:
					#	ref=''.join([ref,*[str(v)+k for k,v in b.items()] ])
					#else:
					#	ref=''.join([ref,"(",*[str(v)+k for k,v in b.items()], ")" ])
					
					if b and max(b.values())>=cov_min-1:
						#cell_cov[bc_strand[0]][0] += sum(b.values())
						#cell_cov[bc_strand[0]][1] += 1
						if (len(b)==2 and min(b.values())==1): # or (print_var and pos==print_var[1] and chrom==print_var[0])
							print(bc_strand[0][5:], chrom, bc_strand[1], pos, *[k for k,v in sorted(b.items(), key=lambda item: item[1], reverse=True)], ''.join([str(v)+k for k,v in b.items()]), ''.join(pcr_ers)) # ref.rstrip("(")[:-1]
							cell_cov[bc_strand[0]][0] += sum(b.values())
							cell_cov[bc_strand[0]][1] += 1
							cell_cov[bc_strand[0]][2]+=1
			#out[chrom][bc_strand][start] = ref

pattern = re.compile(r'(\d+)([MIDNS])') 
pattern_pcr_vars = re.compile(r'(\d+)([ACTG])') 

max_chr_len = 195154279
l='ACGT'
umi2code_d = dict(zip([''.join([a,b,c,d]) for a in l for b in l for c in l for d in l], range(256)))
umi2code_d.update(dict(zip([''.join([a,b]) for a in l for b in l], range(257, 257+16)))) #6)))) #
code2umi = {v: k for k, v in umi2code_d.items()}

out = dict()
cell_cov=dict()
chrom_old=''

for line in sys.stdin:
	line = line.strip().split() 
	if line[-2]=="CB:Z:-" or line[-1]=="UB:Z:-" or line[4]!="255" or line[2]=="chrM" or int(line[1])>=256:
		continue
	chrom = line[2]
	bc = line[-2] ; umi= line[-1] ; strand = int(line[1][0])
	umi = umi2code_d[umi[5:9]]*1000000 + umi2code_d[umi[9:13]]*1000 + umi2code_d[umi[13:17]] + max_chr_len # max size of chr in mice, so it would not intersect with pos after umi_consensus in out
	#if umi == 414217407:
	#	print("+1")
	if chrom != chrom_old:
		if chrom_old and chrom_old in out:
			across_chrom(chrom_old)
			del out[chrom_old]
		chrom_old = chrom
		unpack_startsends(chrom) 
		#print("unpacked", out[chrom_old][('TGTTCCGAGCCTGACC',1)][414217407])

	if (chrom not in out) or ((bc,strand) not in out[chrom]) or (umi not in out[chrom][bc,strand]):
		#if umi == 414217407:
		#	print("not in out")
		continue

	start = int(line[3]); cigar = line[5]
	seq = line[9]

	pileup(start, cigar, seq)
	#if umi == 414217407:
	#	print(out[chrom][(bc,strand)][umi])

with open(address_er, "w") as f:
	for bc in cell_cov:
		if not cell_cov[bc][1]:
			continue
		print(bc[5:], 
			round(cell_cov[bc][2]/cell_cov[bc][1], 8), 
			*cell_cov[bc][::-1],
			file=f)