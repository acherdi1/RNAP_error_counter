#!/usr/bin/env bash

#### USAGE ####
#    This tool calculates transcription error rate 
#    in the given 10x dataset.
#    Stages:
#    0 - download sample (currently under development)
#    1 - align sample to genome
#    2 - subset covered enough positions
#    3 - collect info about substitutions and
#        calculate transcription error rate
#    4 - collect info about error rate genewise

#### PARAMETERS ####
# if you want any of the parameters to be input as a part of terminal command 
# rather than here, e.g. ./main.sh <stage>, replace its value here with $<pos>, 
# where pos is position of the parameter in the command, 
# e.g. ./main.sh <stage> => stage=$1

stage=1

# you can uncomment following variable assignments 
# and variations of other commands linking to them
# in order to automatize assignment for some variables below:

#stage=$1
# dataset=$2
# sample=$3

workdir='/cellfile/datapublic/acherdi1/rnaspeed/'
gitdir="$workdir/RNAP_error_counter/"
out="GSE232273/res/SRR24507709/"
# out="${dataset}/res/${sample}/"

### params for stage 1

# if you have a specific conda environment for running STAR, 
# uncomment following line and adjust it accordingly:
# source /data/public/apapada1/Conda/anaconda/bin/activate star

# path to STAR reference (folder)
genome="/data/public/apapada1/Felix_scRNAseq/starconda" 

# path to the reads from transcripts (if multiple lanes, input all separated by comma):
# please, rewrite it according to the path to the sample of interest:
fastq_read="GSE232273/src/SRR24507709_L1_2.fastq.gz,GSE232273/src/SRR24507709_L2_2.fastq.gz"
#fastq_read="${dataset}/src/${sample}_2.fastq.gz"

# path to the reads from BC+UMIs:
fastq_umi="GSE232273/src/SRR24507709_L2_1.fastq.gz,GSE232273/src/SRR24507709_L2_1.fastq.gz"
# fastq_umi="${dataset}/src/${sample}_1.fastq.gz"

# CellRanger v3, or v2, or custom set-up(0)? (different whitelists for cell barcodes, different UMI lengths)
# https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html
cellranger=2

if (( $cellranger == 3 ))
then 
	# path to the barcode whitelist 
	# (V3: https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz )
	CBWL="$workdir/CBwhitelist.CRv3.3M-february-2018.txt" 
	# length of BC+UMI reads
	umiread_len=150 
	umi_len=12
	umi_params="--soloCBwhitelist $CBWL --soloBarcodeReadLength $umiread_len --soloUMIlen $umi_len"

else if (( $cellranger == 2 ))
then
	# path to the barcode whitelist 
	# (V2: https://teichlab.github.io/scg_lib_structs/data/737K-august-2016.txt.gz )
	CBWL="$workdir/CBwhitelist.CRv2.737K-august-2016.txt"
	# length of BC+UMI reads
	umiread_len=150 
	umi_len=10
	umi_params="--soloCBwhitelist $CBWL --soloBarcodeReadLength $umiread_len --soloUMIlen $umi_len"

else if (( $cellranger == 0 ))
then
	umi_params="--soloBarcodeReadLength 150 --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8"
fi

threads=16

# if you want to add any extra parameters to STAR, then specify them here:
star_params=""
# if you have issues with fitting into RAM limitations of your workstation,
# check up these pages:
# https://github.com/alexdobin/STAR/issues/1159
# https://github.com/alexdobin/STAR/issues/457
# and fiddle these parameters in $star_params:
# star_params="--outBAMsortingThreadN 1 --limitBAMsortRAM 94100100100 --outBAMsortingBinsN 3000


### params for stage 2
# (works only on BAM files sorted by coordinates, 
# that happens automatically in STAR Solo;
# https://github.com/alexdobin/STAR/issues/1792
# processing of unsorted BAM files takes up too much time and RAM,
# at the moment we are improving it)

# path to the BAM file (output of STAR Solo)
bam="${out}/Aligned.sortedByCoord.out.bam"

# path to filtered cell barcodes
bc_filtd="${out}/Solo.out/GeneFull/filtered/barcodes.tsv"

# min amount of duplicates of the same UMI read per position (to account for PCR errors)
th_dup_s=5

# min amount of UMI reads per position (to exclude allele differences)
th_cov_s=5


### params for stage 3
# (works only on BAM files sorted by coordinates, 
# that happens automatically in STAR Solo;
# https://github.com/alexdobin/STAR/issues/1792
# processing of unsorted BAM files takes up too much time and RAM,
# at the moment we are improving it)

# path to the BAM file (output of STAR Solo)
bam="${out}/Aligned.sortedByCoord.out.bam"

# path to subsetting output
pos_filtd="${out}/covd_${th_dup}${th_cov}.txt"

# min amount of duplicates of the same UMI read per position (to account for PCR errors)
th_dup_e=5

# min amount of UMI reads per position (to exclude allele differences)
th_cov_e=5

# allowed amount of mismatches between UMI reads at the position 
# to consider it as transcription error:
mm_allowed=1
# still, some part could be due to allele differences,
# although most of the samples are done on homozygous organisms.
# also, some positions could be more prone to transcription errors than other,
# then there could be more than one RNA carrying the error in the cell,
# and these positions would be excluded.

# output paths
subs_out="${out}/subs_${th_dup}${th_cov}.txt" 
er_out="${out}/er_${th_dup}${th_cov}.txt"

### params for stage 4
subs_out="${out}/subs_${th_dup}${th_cov}.txt"

# path to gene info
gene_info="${workdir}/genes.gex_mm10_2020_A"
# cat genes.gtf |  # /data/public/apapada1/Felix_scRNAseq/refdata-gex-mm10-2020-A/genes/genes.gtf
# awk '$3~/gene/{print $1, $7, $4, $5, substr($16, 2, length($16) - 3)}' > 
# genes.gex_mm10_2020_A
# =>
# chr1 - 3205901 3671498 ENSMUSG00000051951 Xkr4
# chr1 + 3466587 3513553 ENSMUSG00000089699 Gm1992



#### CODE ####

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

if (( $stage < 2 ))
then 
	echo "Running STAR:"
	echo """
	genome=$genome
	CBWL=$CBWL
	umi_len=$umi_len
	threads=$threads
	fastq_read=$fastq_read
	fastq_umi=$fastq_umi
	star_params=$star_params
	""" 

	date
	/usr/bin/time -v STAR --runThreadN $threads --genomeDir $genome \
	--readFilesIn $fastq_read $fastq_umi --soloType Droplet \
	--soloUMIfiltering MultiGeneUMI $umi_params \
	--soloCBmatchWLtype 1MM_multi_pseudocounts --outSAMtype BAM \
	Unsorted SortedByCoordinate --outFileNamePrefix $out \
	--soloFeatures GeneFull --outSAMattributes CB UB GX GN \
	--readFilesCommand zcat $star_params #&> $starlog
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
	/usr/bin/time -v python3 $gitdir/subset_pos.py $th_dup_s $th_cov_s $bc_filtd $pos_filtd
	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Covered positions subsetted"
fi


if (( $stage < 4 ))
then 
	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Calculating RNA Pol 2 error rate:"
	samtools view $bam |
	/usr/bin/time -v python3 $gitdir/er.v1.py \
	$th_dup_e $th_cov_e $mm_allowed $pos_filtd $subs_out $er_out
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
	cat $subs_out | 
	awk 'BEGIN{s[0]="+";s[1]="-"} {print $1, $3, s[$2], $4, $5, $6}' > $out/subs3_${th_dup}${th_cov}.txt
	# sed 's/ 1 / - /g' | sed 's/ 0 / + /g'  > $out/subs3.txt

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "Adding gene info to subs3.txt "
	python3 $gitdir/subs2gene.py $out/subs3_${th_dup}${th_cov}.txt $gene_info \
	> $out/subs_${th_dup}${th_cov}.wg.txt

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "count_vars_per_pos "
	python3 $gitdir/count_vars_per_pos.py $out/subs_${th_dup}${th_cov}.wg.txt \
	> $out/subs_${th_dup}${th_cov}.wg.grouped.txt

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "gene_cov "
	samtools view $bam |
	python3 $gitdir/gene_cov.py \
	$out/subs_${th_dup}${th_cov}.wg.txt $gene_info \
	$out/gene_cov_${th_dup}${th_cov}.txt $out/gene_absent_${th_dup}${th_cov}.txt

	date +"%a %d.%b %T" | tr '\n' ' ' 
	echo "sort "
	cat $out/gene_cov_${th_dup}${th_cov}.txt | sort -rnk 5 > $out/gene_cov_${th_dup}${th_cov}.sorted.txt
fi
echo
echo "DONE!"
date
