#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

PROBEFILE1=$1
PROBEFILE2=$2
PROBEFILE1NAME=$(echo $1 | cut -d . -f 1)
PROBEFILE2NAME=$(echo $2 | cut -d . -f 1)

if [ -d '4probe-reduction/' ]; then

	printf "Starting blast of $PROBEFILE2 against $PROBEFILE1\n"

else 

	mkdir 4probe-reduction

	printf "Starting blast of $PROBEFILE2 against $PROBEFILE1\n"

fi


if [ -f  '/4probe-reduction/'$PROBEFILE2NAME'-'$PROBEFILE1NAME'.blastn' ]; then

	printf "$PROBEFILE2 already BLAST'd against $PROBEFILE1\n"

elif [ -f '0data/probes/'$PROBEFILE1'.nin' ]; then

	cd 0data/probes

	blastn -query $PROBEFILE2 \
	-db $PROBEFILE1 \
	-num_threads 3 \
	-max_target_seqs 100000 \
	-perc_identity 0.90 \
	-outfmt "6 qseqid qstart qend sseqid sstart send gaps mismatch pident evalue bitscore" \
		> '../../4probe-reduction/'$PROBEFILE2NAME'-'$PROBEFILE1NAME'.blastn';	
	cd ../../

	printf "creating list of $PROBEFILE2 loci to remove\n"

	awk '$7 == "0" {print $1}' \
	'4probe-reduction/'$PROBEFILE2NAME'-'$PROBEFILE1NAME'.blastn' | \
	cut -d _ -f 1 | \
	sort -u \
	> '4probe-reduction/'$PROBEFILE2NAME'.to-remove'

else

	printf "Creating $PROBEFILE1NAME database\n"
	
	cd 0data/probes

	makeblastdb -in $PROBEFILE1 -dbtype nucl

	blastn -query $PROBEFILE2 \
	-db $PROBEFILE1 \
	-num_threads 3 \
	-max_target_seqs 100000 \
	-perc_identity 0.90 \
	-outfmt "6 qseqid qstart qend sseqid sstart send gaps mismatch pident evalue bitscore" \
		> '../../4probe-reduction/'$PROBEFILE2NAME'-'$PROBEFILE1NAME'.blastn';	
	cd ../../

	printf "creating list of $PROBEFILE2 loci to remove\n"

	awk '$7 == "0" {print $1}' \
	'4probe-reduction/'$PROBEFILE2NAME'-'$PROBEFILE1NAME'.blastn' | \
	cut -d "p" -f 1 | \
	sort -u \
	> '4probe-reduction/'$PROBEFILE2'.to-remove'

fi
