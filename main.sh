#!/usr/bin/env bash

## CURRENT ISSUES AND TASKS
# NOT ONLY 0 and 16 in $2 => exclude >256, dissociate from chrom variable
# introduce download

#### PARAMETERS ####
dataset=$1 #"GSE232273"
sample=$2 #"SRR24507709"
stage=$3
#0=download datasets and further
#1=STAR and further
#2=subsetting and further
#3=errorrate and further
#4=geneinfo 

r_read="${dataset}/src/${sample}_${4}.fastq.gz"
r_umi="${dataset}/src/${sample}_${5}.fastq.gz"
star_params="$6"
## CellRanger v3 or v2? if not v3, change:
CBWL="CBwhitelist.CRv3.3M-february-2018.txt"
#CBWL="CBWL.CRv2.737K-august-2016.txt"
## which UMI length? if not 12, change:
umi_len=12 #10
threads=16

### not to change 
workdir='/cellfile/datapublic/acherdi1/rnaspeed/'
gitdir="$workdir/RNAP_error_counter/"
out="${dataset}/res/${sample}/"
bam="${out}/Aligned.sortedByCoord.out.bam"
bc_filtd="${out}/Solo.out/GeneFull/filtered/barcodes.tsv"
th_dup=5
th_dup_same=$(($th_dup-1)) #min nr of dups of consensus
th_cov=5 # min nr of reads covering the position
vars_exp=2 # expected nr of variants at the position to consider it as RNAP error
# more - could be biological meaning, not error
gene_info='genes.gex_mm10_2020_A'

star_conda="/data/public/apapada1/Conda/anaconda/bin/activate star"
genome="/data/public/apapada1/Felix_scRNAseq/starconda"
#starlog="${out}/STAR.log"

#### ####

date

echo """
workdir=$workdir
gitdir=$gitdir
dataset=$dataset
sample=$sample
stage=$stage
out=$out

"""

cd $workdir
mkdir -p $out


if (( $stage < 2 )) #star.sh 
then 
	echo "Running STAR:"

	echo """
	star_conda=$star_conda
	genome=$genome
	CBWL=$CBWL
	umi_len=$umi_len
	threads=$threads
	r_read=$r_read
	r_umi=$r_umi
	star_params=$star_params
	""" # starlog=$starlog

	source $star_conda

	date
	/usr/bin/time -v STAR  --runThreadN $threads --genomeDir $genome \
	--readFilesIn $r_read $r_umi --soloType Droplet \
	--soloCBwhitelist $CBWL --soloUMIfiltering MultiGeneUMI \
	--soloCBmatchWLtype 1MM_multi_pseudocounts --outSAMtype BAM \
	Unsorted SortedByCoordinate --outFileNamePrefix $out \
	--soloFeatures GeneFull --outSAMattributes CB UB GX GN \
	--readFilesCommand zcat --soloUMIlen $umi_len $star_params #&> $starlog
	date; echo "STAR finished."
fi
#exit()

if (( $stage < 3 ))
then
	echo "Reads subsetting:"
	echo """
	bam=$bam
	bc_filtd=$bc_filtd
	th_dup=$th_dup
	th_cov=$th_cov
	gene_info=$gene_info
	
	""" 
	
	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Subsetting covered positions (min. $th_cov UMIs per pos, min. $th_dup duplicates per UMI)"
	samtools view $bam |
	awk '$5~/255/ && $2<256 && $(NF-1)!~/CB:Z:-/ && $NF!~/UB:Z:-/' |
	/usr/bin/time -v python3 $gitdir/subset_pos.py $th_dup $th_cov $bc_filtd ${out}/covd.txt
	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Covered positions subsetted"

	#echo "Line numbers extraction:"; date
	#samtools view $bam |
	#/usr/bin/time -v python3 $gitdir/line_nrs.py ${out}/bcumi_dupd.txt > ${out}/line_nrs.txt
	#date; echo "Line numbers are extracted."

	#echo "Lines extraction:"; date
	#samtools view $bam |
	#/usr/bin/time -v awk 'NR==FNR {
	#        linesToPrint[$0]; next
	#} FNR in linesToPrint {
	#        if($(NF-2)~/^G[NX]:Z:/){gene=substr($(NF-2),6)};
	#        print substr($(NF-1),6), $3"_"substr($2,1,1), 
	#        substr($NF,6), $4, $6, $10, gene
	#}' ${out}/line_nrs.txt - > ${out}/reads.unsorted.txt
	#/usr/bin/time -v sort ${out}/reads.unsorted.txt > ${out}/reads.txt && rm ${out}/reads.unsorted.txt
	#date; echo "Lines are extracted."

	#echo "Extraction of UMIs covering at least 1NT equal or more than $th_cov times:"; date
	#/usr/bin/time -v python3 $gitdir/bcumi_covd.py \
	#$th_dup_same $th_cov $out/reads.txt $out/bcumi_covd.txt
	#date; echo "Barcodes and UMIs are extracted."
fi


if (( $stage < 4 ))
then 
	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Calculating RNA Pol 2 error rate:"
	samtools view $bam |
	/usr/bin/time -v python3 $gitdir/er.v2.py \
	$th_dup $th_cov $vars_exp $out/covd.txt \
	$out/subs.txt $out/er.txt
	date +"%a %d.%b %T"| tr '\n' ' ' 
	echo "Error rate is calculated."
fi

if (( $stage < 5 ))
then 
	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Introducing gene information:"

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Editing subs.txt (adding strand info) "
	# rewrite er.v2.py and exclude this first chunk
	cat $out/subs.txt | 
	awk 'BEGIN{s[0]="+";s[1]="-"} {print $1, $3, s[$2], $4, $5, $6}' > $out/subs3.txt
	# sed 's/ 1 / - /g' | sed 's/ 0 / + /g'  > $out/subs3.txt
	#date; echo "Edited (=>subs3.txt)."

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Adding gene info to subs3.txt "
	python3 $gitdir/subs2gene.py $out/subs3.txt $gene_info \
	> $out/subs.wg.txt
	# date; echo "Added(subs.wg.txt)."

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "count_vars_per_pos "
	python3 $gitdir/count_vars_per_pos.py $out/subs.wg.txt \
	> $out/subs.wg.grouped.txt
	#date | tr '\n' ' ' ; echo "Counted."

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "gene_cov "
	samtools view $bam |
	python3 $gitdir/gene_cov.py \
	$out/subs.wg.txt $gene_info \
	$out/gene_cov.txt $out/gene_absent.txt

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "sort "
	cat $out/gene_cov.txt | sort -rnk 5 > $out/gene_cov.sorted.txt
fi
echo
echo "DONE!"
date
