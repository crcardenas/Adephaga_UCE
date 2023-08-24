#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

GEN=$(echo $1 | cut -d . -f 1)
PROBE=$(echo $2 | cut -d . -f 1)

# index genome and map probes
if [ -f '6group-by/'$PROBE'-'$GEN'.sam' ]; then

	printf "$PROBE already mapped to -$GEN \n";

elif [ -f '0data/genomes/'$GEN'.fna.pac' ]; then

	printf "Mapping $PROBE to $GEN \n";

	bwa mem '0data/genomes/'$GEN'.fna' '5join-probes/'$PROBE'.fasta' \
		> '6group-by/'$PROBE'-'$GEN'.sam';

else

	printf "Indexing $GEN \n";
	bwa index '0data/genomes/'$GEN'.fna';
	
	printf "Mapping $PROBE to $GEN \n";
	bwa mem '0data/genomes/'$GEN'.fna' '5join-probes/'$PROBE'.fasta' \
		> '6group-by/'$PROBE'-'$GEN'.sam';
		
	printf "finished bwa indexing and mapping\noutput $PROBE-$GEN.sam in 6group-by directory\n";

fi

# convert sam to bam
if [ -f '6groupby/'$PROBE'-'$GEN'.mem.sorted.mapped.bed' ]; then

	printf "BED file already created for $PROBE-$GEN\n"

else

	printf "converting to bam file\n"

	samtools view -h -b -S '6group-by/'$PROBE'-'$GEN'.sam' \
		> '6group-by/'$PROBE'-'$GEN'.mem.bam'

	printf "sorting reads\n"
	
	# sort bam file
	samtools sort '6group-by/'${PROBE}'-'${GEN}'.mem.bam' \
		-o '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.bam'

	printf "extracting $PROBE probes that map to $GEN\n"
	
	# extract sequences that map against the genome
	samtools view -b -F 4 \
		'6group-by/'${PROBE}'-'${GEN}'.mem.sorted.bam' \
		> '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bam'

	printf "creating BED file\n"
	
	# create bed file
	bedtools bamtobed -i '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bam' \
		> '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bed'

	printf "output $PROBE-$GEN.mem.sorted.mapped.bed in 6group-by/ directory\n"

fi

echo -e "Find intersect of $PROBE and $GEN gene features \n"
bedtools intersect -a '1genefeatures/'${GEN}'.sorted.genes.gff' \
    -b '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bed' \
    -names ${PROBE} \
    -wa -wb > '6group-by/'${PROBE}-${GEN}.gene.intersect

#here we group the GFF feature info and sumarize by UCE loci
# -g group by column or ranges (1-4; 1,2,3,4)
# -c summarize by this column
printf "group $PROBE mapped to $GEN \n"
bedtools groupby -i '6group-by/'${PROBE}-${GEN}.gene.intersect\
	 -g 9 -c 13 -o collapse > '6group-by/'${PROBE}'-'${GEN}'.grpby_gene'

# create output files

cat '6group-by/'${PROBE}'-'${GEN}'.grpby_gene' | tr -d " " > 6group-by/tmp
while read TMP; do
	GENE=$(echo $TMP | cut -d";" -f 1 | cut -d: -f 2)
    UCELIST=$(echo $TMP | awk '{print $2}')
    echo -e "$GENE \t $UCELIST"
done < '6group-by/tmp' | tr -d " " > '6group-by/tmp2'

GENOME=$(echo $GEN | cut -d"_" -f1,2)
awk 'BEGIN{FS=OFS="\t"} 																		# set intput and output as tab delimited
	{split($2, probes, ","); 																	# split the second feild of each line into an array called probes, and as a comma delimited file
	for (i in probes) {split(probes[i], loci, "_"); 											# loop through each element in the proe array, and split it according to the underscore
	loci_array[loci[1]] = 1} delete probes; 													# while in the for loop save the first part of the loci as associative array called loci_array
	for (j in loci_array) {loci_list = loci_list ? loci_list "," j : j}; 						# loop through each key in the loci array and append a string called "loci list". The string is built up using the ? : in the loop condition: if loci_list is not empty (if it does not contain a previous loci) the current loci is appended with a comma; otherwise the current loci is added to the string.
	delete loci_array; 																			# delete the loci array to free up memory for the next line
	print $1, loci_list; loci_list=""}' 6group-by/tmp2 > '6group-by/'${PROBE}'-'${GEN}'.short.grpby_gene'		# print the first field ($1) followed by the loci_list string built earlier and then reset loci_list to start with the next line.

printf "finished grouping $PROBE by $GEN genes \n"
printf "\t Detailed grouped by file is ${PROBE}-${GEN}.grpby_gene\n"
printf "\t Succinct grouped by gene file is ${PROBE}-${GEN}.short.grpby_gene\n"
