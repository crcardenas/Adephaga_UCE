#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

# turn off N probes (e.g., PROBE1:PROBEN) depending on probesets used for characterization

GEN=$1
PROBE1=$2
PROBE2=$3
PROBE3=$4

printf "running index script\n"

# index probes against genomes
bash index.sh $GEN $PROBE1;
bash index.sh $GEN $PROBE2;
bash index.sh $GEN $PROBE3;

printf "running sam2bed script\n"

# convert sam file to bam file
bash sam2bed.sh $GEN $PROBE1;
bash sam2bed.sh $GEN $PROBE2;
bash sam2bed.sh $GEN $PROBE3;

printf "ready to run intersect.sh\n"