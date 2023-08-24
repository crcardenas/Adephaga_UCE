#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

GEN=$(echo $1 | cut -d . -f 1)
PROBE=$(echo $2 | cut -d . -f 1)

# use BWA to map using mem
# check for indexed genome files and if none
# are present index them and then use
if [ -f '2mapping/'$PROBE'-'$GEN'.sam' ]; then
	printf "$PROBE already mapped to -$GEN \n\n";

elif [ -f '0data/genomes/'$GEN'.fna.pac' ]; then  
	printf "Mapping $PROBE to $GEN \n";
	bwa mem '0data/genomes/'$GEN'.fna' '0data/probes/'$PROBE'.fasta' \
		> '2mapping/'$PROBE'-'$GEN'.sam';

else
	printf "Indexing $GEN \n";
	bwa index '0data/genomes/'$GEN'.fna';
	
	printf "Mapping $PROBE to $GEN \n";
	bwa mem '0data/genomes/'$GEN'.fna' '0data/probes/'$PROBE'.fasta' \
		> '2mapping/'$PROBE'-'$GEN'.sam';
	
	printf "Finished bwa indexing and mapping\nOutput $PROBE-$GEN.sam in 2mapping directory \n\n";
fi