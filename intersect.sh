#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

# turn off N probes (e.g., PROBE1:PROBEN) depending on probesets used for characterization

# user bedfile input
BED1=$2
BED2=$3
BED3=$4
# create variable for naming probes and genomes
GEN=$(echo $1 | cut -d . -f 1)
PROBE1=$(echo $2 | cut -d . -f 1 | cut -d - -f 1)
PROBE2=$(echo $3 | cut -d . -f 1 | cut -d - -f 1)
PROBE3=$(echo $4 | cut -d . -f 1 | cut -d - -f 1)

printf "getting intersect of all probes with $GEN \n"
# all probes intersect
bedtools intersect \
	-a  '1genefeatures/'$GEN'_simple.sorted.gff' \
	-b '2mapping/'$BED1 \
	'2mapping/'$BED2 \
	'2mapping/'$BED3 \
	-names $PROBE1 $PROBE2 $PROBE3 \
	-wa -wb > '3intersect/'$GEN'-PROBES.intersect';
# turn off N -b inputs depending on the total number of probes you use
#	-names $PROBE2 $PROBE3 \
# summary
bedtools summary -i '3intersect/'$GEN'-PROBES.intersect' \
	-g '0data/genomes/'$GEN'.genomefile' \
	> '3intersect/'$GEN'-PROBES-intersect-summary.txt'

printf "getting which probes do not intersect with $PROBE1\n"
# get all probes that do not intersect with probe1
bedtools intersect \
	-a '2mapping/'$BED1 \
	-b '2mapping/'$BED2 \
	'2mapping/'$BED3 \
	-names $PROBE2 $PROBE3 \
	-v > '3intersect/'$GEN'-PROBES.nointersect';

# summary
bedtools summary -i '3intersect/'$GEN'-PROBES.nointersect' \
	-g '0data/genomes/'$GEN'.genomefile' \
	> '3intersect/'$GEN'-nointersect-summary.txt'

printf "getting which probes do not intersect with $PROBE2\n"
# get all probes that do not intersect with probe1
bedtools intersect \
	-a '2mapping/'$BED2 \
	-b '2mapping/'$BED1 \
	'2mapping/'$BED3 \
	-names $PROBE1 $PROBE3 \
	-v > '3intersect/'$GEN'-'$PROBE2'.nointersect';

# summary
bedtools summary -i '3intersect/'$GEN'-'$PROBE2'.nointersect' \
	-g '0data/genomes/'$GEN'.genomefile' \
	> '3intersect/'$GEN'-'$PROBE2'-nointersect-summary.txt'

printf "getting which probes do not intersect with $PROBE3\n"
# get all probes that do not intersect with probe3
bedtools intersect \
	-a '2mapping/'$BED3 \
	-b '2mapping/'$BED2 \
	'2mapping/'$BED1 \
	-names $PROBE2 $PROBE1 \
	-v > '3intersect/'$GEN'-'$PROBE3'.nointersect';

# summary
bedtools summary -i '3intersect/'$GEN'-'$PROBE3'.nointersect' \
	-g '0data/genomes/'$GEN'.genomefile' \
	> '3intersect/'$GEN'-'$PROBE3'-nointersect-summary.txt'

printf "intersect of probes complete\n"