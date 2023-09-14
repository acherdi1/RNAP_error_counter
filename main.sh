#!/usr/bin/env bash

## CURRENT ISSUES AND TASKS
# NOT ONLY 0 and 16 in $2 => exclude >256, dissociate from chrom variable
# MOVE er.p1.fin.py into bcumis_mt_th.py
# introduce download
# extend th_dup and th_cov into scripts
# rewrite scripts so i could see all the names here

#### PARAMETERS ####
dataset=$1 #"GSE232273"
sample=$2 #"SRR24507709"
stage=$3
#0=download datasets and further
#1=STAR and further
#2=subsetting and further
#3=errorrate and further
#4=geneinfo 

### not to change 
workdir='/cellfile/datapublic/acherdi1/rnaspeed/'
gitdir="$workdir/RNAP_error_counter/"
out="${dataset}/res/${sample}/"
bam="${out}/Aligned.sortedByCoord.out.bam"
bc_filtd="${out}/Solo.out/GeneFull/filtered/barcodes.tsv"
th_dup=5
th_cov=5
gene_info='genes.gex_mm10_2020_A'

#### STAR settings, change only if run STAR ####

## CellRanger v3 or v2? if not v3, change:
CBWL="CBwhitelist.CRv3.3M-february-2018.txt" #"CBWL.CRv2.737K-august-2016.txt"

## which UMI length? if not 12, change:
umi_len=12 #10

r_read="${dataset}/src/${sample}_4.fastq.gz"
r_umi="${dataset}/src/${sample}_3.fastq.gz"
threads=16

### not to change
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

	""" # starlog=$starlog

	source $star_conda

	date
	/usr/bin/time -v STAR  --runThreadN $threads --genomeDir $genome \
	--readFilesIn $r_read $r_umi --soloType Droplet \
	--soloCBwhitelist $CBWL --soloUMIfiltering MultiGeneUMI \
	--soloCBmatchWLtype 1MM_multi_pseudocounts --outSAMtype BAM \
	SortedByCoordinate --outFileNamePrefix $out \
	--soloFeatures GeneFull --outSAMattributes CB UB GX GN \
	--readFilesCommand zcat --soloUMIlen $umi_len  &> $starlog
	date; echo "STAR finished."
fi

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
	
	echo "Extraction of cell barcodes and UMIs present in BAM file more than $th_dup times:"; date
	samtools view $bam |
	/usr/bin/time -v python3 $gitdir/bcumi_mt_th.py $th_dup $bc_filtd > ${out}/bcumis_mt_th.txt
	date; echo "Barcodes and UMIs are extracted."

	echo "Line numbers extraction:"; date
	samtools view $bam |
	/usr/bin/time -v python3 $gitdir/line_nrs.py ${out}/bcumis_mt_th.txt > ${out}/line_nrs.txt
	date; echo "Line numbers are extracted."

	echo "Lines extraction:"; date
	samtools view $bam |
	/usr/bin/time -v awk 'NR==FNR {
	        linesToPrint[$0]; next
	} FNR in linesToPrint {
	        if($(NF-2)~/^G[NX]:Z:/){gene=substr($(NF-2),6)};
	        print substr($(NF-1),6), $3"_"substr($2,1,1), 
	        substr($NF,6), $4, $6, $10, gene
	}' ${out}/line_nrs.txt - | /usr/bin/time -v sort > ${out}/reads.txt
	date; echo "Lines are extracted."

	echo "Extraction of UMIs covering at least 1NT equal or more than $th_dup times:"; date
	/usr/bin/time -v python3 $gitdir/er.p1.fin.py $out/
	date; echo "Barcodes and UMIs are extracted."
fi


if (( $stage < 4 ))
then 
	echo "Calculating RNA Pol 2 error rate:"; date
	/usr/bin/time -v python3 $gitdir/er.p2.fin.py $out/
	date; echo "Error rate is calculated."
fi

if (( $stage < 5 ))
then 
	echo "Introducing gene information:"

	echo "Editing subs.txt (adding strand info)"; date
	cat $out/subs.txt | 
	sed 's/_1/ -/g' | sed 's/_0/ +/g'  > $out/subs3.txt
	date; echo "Edited (subs3.txt)."

	echo "Adding gene info to subs3.txt"; date
	python3 $gitdir/subs2gene.py $out/subs3.txt $gene_info \
	> $out/subs3.wg.txt
	date; echo "Added(subs3.wg.txt)."

	echo -ne "count_vars_per_pos "; date
	python3 $gitdir/count_vars_per_pos.py $out/subs3.wg.txt \
	> $out/subs3.wg.grouped.txt

	echo -ne "gene_cov "; date
	python3 $gitdir/gene_cov.py $out/subs3.wg.txt \
	$out/reads.txt $gene_info \
	$out/gene_cov.txt $out/gene_absent.txt

	echo -ne "sort "; date
	cat $out/gene_cov.txt | sort -rnk 5 > $out/gene_cov.sorted.txt
fi
