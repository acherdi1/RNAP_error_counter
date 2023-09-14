import re
import sys
import numpy as np
from collections import defaultdict

dup_min = 4
cov_min = 5
vars_exp = 2
work_dir = sys.argv[1]
input_file = '/'.join([work_dir, "reads.txt"])
bcumis_file='/'.join([work_dir, "bcumis.me5dup_me5cov.txt"])
address_subs='/'.join([work_dir, "subs.txt"])
address_er='/'.join([work_dir, "er.txt"])

with open(address_subs, 'w') as f:
	pass
with open(address_er, 'w') as f:
	pass

# maybe there should be min nr of covered nts to consider a cell covered?
# now the whole statistics of errors is counted based on
# overall nr of cells as those cells where at least 1 nt is covered
# is it a proper way?


# вместо словаря на атцг сделать список нт,колво, и сверять с нт, 
# если нт другой то сразу дропать, если этот то плюсовать к колву
# тогда без внутр словарей мб лучше будет чиститься память?

# НО ГЛАВНОЕ СЕЙЧАС УСКОРЕНИЕ!!

#разделить сикв на отрезки по 4 => 
#=> словарь AAAA:0, AAAC:1...TTTT:255 =>
#сравнивать не нт а отрезки нт, и 
#хранить только отрезки нт
#работало бы, не будь разрывов. 
#С разрывами надо как-то умгеен

# if on umi level i would have another dict, then i could have introduced consensus sequences to out dict already, 
# in a right nestness, and by that i would not need to do umipos2posumi

# np.zeros in pileup to another type?

# per cell:
# motifs=neighbouring nts if they are covered 5 or more times
# in R script, minus strands are rev, but redirected in bam!!!


def pileup(start, cigar, seq_len):
	pos = int(start)
	pos_in = 0
	coords = np.zeros(seq_len, dtype=np.uint32)
	cigar_ops = pattern.findall(cigar) 
	if cigar_ops[0][1] == 'S':
		pos_in += int(cigar_ops[0][0])
		cigar_ops = cigar_ops[1:]
	for length, op in cigar_ops:
		length = int(length)
		if op == 'M':
			coords[pos_in:(pos_in+length)] = range(pos, pos + length)
			pos_in += length
		elif op == 'I':
			pos_in += length
			continue
		pos += length
	return(coords)

def leave_only_max(d):
	dmax=[0,0]
	for i in d:
		if d[i]>dmax[1]:
			dmax=[i,d[i]]
	return(dmax)

def proc1(bc, chrom, umi, start, cigar, seq, seq_len):
	coords = pileup(start, cigar, seq_len)
	seq = [nt2int[i] for i in list(seq)]
	for i in range(seq_len):
		out[bc][chrom][umi][coords[i]][seq[i]] += 1
	return()

def umi_consensus(bc, chrom, umi): #, dup_min
	for pos in list(out[bc][chrom][umi].keys()):
		nt, cov = leave_only_max(out[bc][chrom][umi][pos])
		if cov < dup_min or not nt:
			del out[bc][chrom][umi][pos]
		else:
			out[bc][chrom][umi][pos] = nt

	if not out[bc][chrom][umi]:
		del out[bc][chrom][umi]

	return() 

def umipos2posumi(bc, chrom):
	for umi in list(out[bc][chrom].keys()):
		for pos in list(out[bc][chrom][umi].keys()):
			out[bc][chrom][pos][umi] = out[bc][chrom][umi][pos]
			del out[bc][chrom][umi][pos]
		del out[bc][chrom][umi]
	
	if out[bc][chrom][0]:
		del out[bc][chrom][0]

	return()

def variance(bc, chrom): 
	len_cov_chr = 0
	len_cov_err_chr = 0
	for pos in list(out[bc][chrom].keys()):
		if len(out[bc][chrom][pos]) < cov_min:
			del out[bc][chrom][pos]
		else:
			len_cov_chr += 1
			d = dict() 
			for umi in list(out[bc][chrom][pos].keys()):
				try:
					d[out[bc][chrom][pos][umi]] += 1
				except:
					d[out[bc][chrom][pos][umi]] = 1
			variants = len(d)
			if variants != vars_exp or min(d.values())>1:
				del out[bc][chrom][pos]
				# 1) remove those where mismatch is repeated more than once on the same position!
				# so when more than one umi contains mismatch on the same position
				# bc pol2 error hardly could happen twice at the same position
				# meaning that mismatch probably happened due to some other reasons
				# 2) remove if there are more than 2 variants per min position
				# same reasoning
			else:
				len_cov_err_chr += 1
				out[bc][chrom][pos] = sorted(d, key=d.get, reverse=True) #d
				out[bc][chrom][pos] = [int2nt[i] for i in out[bc][chrom][pos]]
				# mismatch from which to which nt
				#[1] " AC" " TA"

	return(len_cov_chr, len_cov_err_chr)

def record_substitutions(address, bc):
	with open(address, "a") as f:
		for chrom in list(out[bc].keys()):
			for coord in list(out[bc][chrom].keys()):
				print(bc, chrom, coord, *out[bc][chrom][coord], file = f)

def record_errorrate(address, bc, er, len_er, len_covd):
	with open(address, "a") as f:
		print(bc, er, len_er, len_covd, file = f)

pattern = re.compile(r'(\d+)([MIDNS])') 
nt2int = {"N":0, "A":1,"C":2,"T":3,"G":4}
int2nt = {v: k for k, v in nt2int.items()}
default_genome = lambda: defaultdict(lambda: defaultdict(
	lambda: defaultdict(lambda: 
		{1: 0, 2: 0, 3: 0, 4: 0, 0: 0})))
out = defaultdict(default_genome)

def load_bcumis(bcumis_file):
	bcumis = list()
	with open(bcumis_file) as f:
		for line in f: 
			bcumis.append(line.strip())
	bcumis.append('bcchromumi') # so not out of range
	return bcumis

def get_seq_len(input_file): 
	with open(input_file) as f:
		seq = f.readline().strip().split()[5]
	seq_len = len(seq)
	return seq_len

def main():
	len_cov = 0; len_cov_err = 0
	bcumis = load_bcumis(bcumis_file)

	seq_len = get_seq_len(input_file)
	old_site = False
	bcumis_pos = 0
	check = False ; bc, chrom, umi = '','',''

	with open(input_file) as f:
		for line in f: 
			#print(line.strip()); print(bcumis[bcumis_pos])

			if not line.startswith(bcumis[bcumis_pos]):
				if check:
					bcumis_pos += 1
					if not line.startswith(bcumis[bcumis_pos]):
						check = False
						continue
				else:
					continue
			check = True
			bc, chrom, umi, start, cigar, seq = line.strip().split()[0:6]
			site = (bc, chrom, umi)

			if site != old_site and old_site: 
				_bc, _chrom, _umi = old_site 
				umi_consensus(_bc, _chrom, _umi) #, dup_min

				if chrom != _chrom or bc != _bc:
					umipos2posumi(_bc, _chrom)
					len_cov_chr, len_cov_err_chr = variance(_bc, _chrom) #, cov_min, vars_exp
					len_cov += len_cov_chr
					len_cov_err += len_cov_err_chr

				if bc != _bc:
					if len_cov>0:
						errorRate=len_cov_err/len_cov
						record_substitutions(address_subs, _bc)
						record_errorrate(address_er, _bc, errorRate, len_cov_err, len_cov)
						len_cov = 0; len_cov_err = 0
					del out[_bc]

			proc1(bc, chrom, umi, start, cigar, seq, seq_len)
			old_site=site

			
		umi_consensus(bc, chrom, umi) #, dup_min
		umipos2posumi(bc, chrom)
		len_cov_chr, len_cov_err_chr = variance(bc, chrom) #, cov_min, vars_exp
		len_cov += len_cov_chr
		len_cov_err += len_cov_err_chr
		if len_cov>0:
			errorRate=len_cov_err/len_cov

			record_substitutions(address_subs, bc)
			record_errorrate(address_er, bc, errorRate, len_cov_err, len_cov)
main()
