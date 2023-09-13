import re
import sys

dup_min = 4
cov_min = 5
work_dir = sys.argv[1]
input_file = '/'.join([work_dir, "reads.txt"])
output_file='/'.join([work_dir, "bcumis.me5dup_me5cov.txt"])
with open(output_file, 'w') as f:
	pass

# subset UMIs with at least 5dups which cover at least 5 times same region

def cigar2pos(start, cigar): 
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
	while startss:
		if starts < ends:
			cov += 1
			if cov >= th and not chunk: #switch_cov5:
				chunk = [starts] #(starts,)
			starts = startss.pop()
		elif starts > ends:
			cov -= 1
			if chunk and cov < th : #switch_cov5
				chunk.append(ends) #+= (ends,)
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
				chunk = (starts,)
			starts=endss[0]
		cov -= 1
		if chunk and cov < th : 
			chunk += (ends,)
			cov5 += (chunk, )
			chunk = False
		ends = endss.pop()
	return(cov5)

pattern = re.compile(r'(\d+)([MIDNS])') 

def main():
	startss = list() #tuple()
	endss = list() #tuple()
	covd = dict()
	with open(input_file) as f:
		old_site = False
		for line in f: #.readlines():
			bc, chrom, umi, start, cigar = line.strip().split()[0:5]
			site = (bc, chrom, umi)
			
			if site != old_site and old_site:
				_bc, _chrom, _umi = old_site
				#if old_site == ('AACTTCTCAAATGCGG', 'chr19', 'AGATCATAGGGG'):
				#	print(sorted(startss))
				#	print(sorted(endss))
				cov5 = covv(startss, endss, th=dup_min)
				#if old_site == ('AACTTCTCAAATGCGG', 'chr19', 'AGATCATAGGGG'):
				#	print(_umi, cov5)
				if cov5:
					covd[_umi] = cov5

				if chrom != _chrom or bc != _bc:
					#if _bc=="AAAGGTAGTCACGCTG" and _chrom=="chr6":
					#	print(covd)
					if len(covd) >= cov_min:
						startss = list()
						endss = list()
						for chunks in covd.values():
							for chunk in chunks:
								startss.append(chunk[0])
								endss.append(chunk[1])

						#startss = sorted(startss, reverse=True)
						#endss = sorted(endss, reverse=True)
						cov5 = covv(startss, endss, th=cov_min)
						#if _bc=="AAAGGTAGTCACGCTG" and _chrom=="chr6":
						#	print(_umi, cov5)
						if cov5:
							with open(output_file, 'a') as f:
								for ummi in covd.keys():
									print(_bc, _chrom,ummi, file=f) #, sep="\t", file=f)

					covd = dict()

				startss = list() #tuple()
				endss = list() #tuple()	

			starts, ends = cigar2pos(start, cigar) # pos
			startss.extend(starts) #+= (*starts,)
			endss.extend(ends) #+= (*ends, )
			old_site = site

main()



#AACTTCTCAAATGCGG	chr19	AGATCATAGGGG		      <
#AACTTCTCAAATGCGG	chr19	AGATGGTGCCCC		      <
#AACTTCTCAAATGCGG	chr19	AGGGTAAGTGTA		      <
#AACTTCTCAAATGCGG	chr19	ATCAACCCTCGT		      <

