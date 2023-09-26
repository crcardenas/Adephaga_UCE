#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

GEN=$(echo $1 | cut -d . -f 1)
PROBE=$(echo $2 | cut -d . -f 1)

# convert sam to bam
if [ -f '2mapping/'$PROBE'-'$GEN'.mem.sorted.mapped.bed' ]; then

	printf "BED file already created for $PROBE-$GEN\n"

else

	printf "converting to bam file\n"

	samtools view -h -b -S '2mapping/'$PROBE-$GEN'.sam' \
		> '2mapping/'$PROBE'-'$GEN'.mem.bam'

	printf "sorting reads\n"
	
	# sort bam file
	samtools sort '2mapping/'$PROBE'-'$GEN'.mem.bam' -o '2mapping/'$PROBE'-'$GEN'.mem.sorted.bam'

	printf "extracting $PROBE probes that map to $GEN\n"
	
	# extract sequences that map against the genome
	samtools view -b -F 4 '2mapping/'$PROBE'-'$GEN'.mem.sorted.bam' \
		> '2mapping/'$PROBE'-'$GEN'.mem.sorted.mapped.bam'

	printf "creating BED file\n"
	
	# create bed file
	bedtools bamtobed -i '2mapping/'$PROBE-$GEN.mem.sorted.mapped.bam \
		> '2mapping/'$PROBE'-'$GEN'.mem.sorted.mapped.bed'

	printf "output $PROBE-$GEN.mem.sorted.mapped.bed in 2mapping/ directory\n"

fi

if [ -f '2mapping/'$PROBE'-'$GEN'.mem.sorted.mapped.bam' ]; then

rm 2mapping/*.bam

printf "Checking for potential loci/probes map to different regions in $PROBE-$GEN.sam\n"

else

printf "Checking for potential loci/probes map to different regions in $PROBE-$GEN.sam\n"

fi

# create simple summary file, create variables and push to variables or new file for: number of UCEs/Loci maped vs UCEs/Loci total, number of probes mapping to different regions.
#number of mapped loci
awk '{if ($3 != "*") print $0}' '2mapping/'$PROBE'-'$GEN'.sam' | \
	grep -v "@" |  \
	awk '{print $1}' | \
	cut -d _ -f 1 | \
	sort -u | \
	grep "." -c \
		> '2mapping/'tmp.uce.count

printf "Total  of $PROBE loci mapped against $GEN\t" && cat 2mapping/tmp.uce.count

# number of multiple mappings
# count multi mapping
cat '2mapping/'$PROBE'-'$GEN'.sam' | 
	grep -v "@" | \
	grep "SA:" -c \
		> 2mapping/tmp.multi-uce.count

printf "Total $PROBE loci with multiple mappings in $GEN\t" && cat 2mapping/tmp.multi-uce.count

# store  names
cat '2mapping/'$PROBE'-'$GEN'.sam' | grep -v "@" | grep "SA:" > 2mapping/tmp.multi

# count UCE/Loci names
awk '{print $1}' 2mapping/tmp.multi | cut -d "p" -f 1 | sort -u | wc -l > 2mapping/tmp.multi-uce.count

# list of uce to remove
awk '{print $1}' 2mapping/tmp.multi | cut -d "p" -f 1 | sort -u > '2mapping/'$PROBE'_to_remove.list'

# store summary
if [ -f 2mapping/summary.txt ]; then

	printf "Writing mapping summary for $PROBE-$GEN\n"

else

	printf "Writing mapping summary for $PROBE-$GEN \n"

	printf "probe-genome\tmapped-loci\ttotal-loci\tNumber-multi-mapped-loci\n" > 2mapping/summary.txt

fi

NAME=$(printf "$PROBE-$GEN")

MAPPEDLOCI=$(cat 2mapping/tmp.uce.count)

TOTALLOCI=$(cat '0data/probes/'$PROBE'.fasta' | grep ">" | cut -d "|" -f 1 | cut -d _ -f1 | sort -u | grep "." -c)

MULTIMAPPING=$(cat 2mapping/tmp.multi-uce.count)

printf "$NAME\t$MAPPEDLOCI\t$TOTALLOCI\t$MULTIMAPPING\n" >> 2mapping/summary.txt

rm 2mapping/tmp.*

# if conditional statement is satisfied then remap with removed loci 

if [ $(cat '2mapping/'$PROBE'-'$GEN'.sam' | grep -v "@" | grep "SA:" -c) == 0 ]; then

	printf "$PROBE does not need loci/probes removed and re-mapped against $GEN \n";


else 

	printf "$PROBE needs re-mapped with removed loci/probes\n"

	# remove loci in probe file based on a list
	awk '{ 
		if ((NR>1) && ($0~/^>/)) { printf("\n%s", $0); }
		else if (NR==1) { printf("%s", $0); } 
		else { printf("\t%s", $0); } 
		}' '0data/probes/'$PROBE'.fasta' | \
	grep -vf '2mapping/'$PROBE'_to_remove.list' - | \
	tr "\t" "\n" \
		> '0data/probes/'$PROBE'_2.fasta';

	printf "Probes removed fom $PROBE\t" && grep ">" '0data/probes/'$PROBE'.fasta' -c

	# use BWA to remap using mem
	# check for indexed genome files and if none are present index them and then map
	if [ -f '2mapping/'$PROBE'-'$GEN'_2.sam' ]; then

		printf "$PROBE already mapped to $GEN \n";

	elif [ -f '0data/genomes/'$GEN'.fna.pac' ]; then  

		printf "mapping new $PROBE set to $GEN \n";

		bwa mem '0data/genomes/'$GEN'.fna' '0data/probes/'$PROBE'_2.fasta' \
			> '2mapping/'$PROBE'-'$GEN'_2.sam';

	else

		printf "indexing $GEN \n";

		bwa index '0data/genomes/'$GEN'.fna';
		
		printf "\n re mapping new $PROBE set to $GEN \n";
		
		bwa mem '0data/genomes/'$GEN'.fna' '0data/probes/'$PROBE'_2.fasta' \
			> '2mapping/'$PROBE'-'$GEN'_2.sam';
		
		printf " new mapped file $PROBE-$GEN_2.sam in 2mapping directory \n";
	
	fi

	# convert sam to bam
	if [ -f '2mapping/'$PROBE'-'$GEN'_2.mem.sorted.mapped.bed' ]; then
		
		printf "BED file already created for new $PROBE-$GEN\n"
	
	else
		
		printf "converting to bam file\n"
		
		samtools view -h -b -S '2mapping/'$PROBE'-'$GEN'_2.sam' \
			> '2mapping/'$PROBE'-'$GEN'_2.mem.bam'
		
		printf "sorting reads\n"

		# sort bam file
		samtools sort '2mapping/'$PROBE'-'$GEN'_2.mem.bam' -o '2mapping/'$PROBE'-'$GEN'_2.mem.sorted.bam'
		
		printf "extracting new $PROBE probes that map to $GEN\n"
		
		# extract sequences that map against the genome
		samtools view -b -F 4 '2mapping/'$PROBE'-'$GEN'_2.mem.sorted.bam' \
			> '2mapping/'$PROBE'-'$GEN'_2.mem.sorted.mapped.bam'
		
		printf "creating BED file\n"
		
		# create bed file
		bedtools bamtobed -i '2mapping/'$PROBE'-'$GEN'_2.mem.sorted.mapped.bam' \
			> '2mapping/'$PROBE'-'$GEN'_2.mem.sorted.mapped.bed'

		printf "output $PROBE-$GEN_2.mem.sorted.mapped.bed in 2mapping/ directory\n"
	
	fi

	if [ -f '2mapping/'$PROBE'-'$GEN'_2.mem.sorted.mapped.bam' ]; then

		rm 2mapping/*.bam
	
		printf "Checking if loci/probes map to different regions in $PROBE-$GEN.sam\n"
	
	else
	
		printf "Checking if loci/probes map to different regions in $PROBE-$GEN.sam\n"
	
	fi

	# rewrite summary into new file
	awk '{if ($3 != "*") print $0}' '2mapping/'$PROBE'-'$GEN'_2.sam' | \
		grep -v "@" |  \
		awk '{print $1}' | \
		cut -d _ -f 1 | \
		sort -u | \
		grep "." -c \
			> '2mapping/'tmp.uce.count

	printf "Total of new $PROBE loci mapped against $GEN\t" && cat 2mapping/tmp.uce.count

	# number of multiple mappings
	# count multi mapping
	cat '2mapping/'$PROBE'-'$GEN'_2.sam' | \
		grep -v "@" | \
		grep "SA:" -c \
			> 2mapping/tmp.multi-uce.count

	printf "Total $PROBE loci with multiple mappings in $GEN\t" && cat 2mapping/tmp.multi-uce.count

	# store  names
	cat '2mapping/'$PROBE'-'$GEN'_2.sam' | \
		grep -v "@" | \
		grep "SA:" \
			> 2mapping/tmp.multi

	# count UCE/Loci names
	awk '{print $1}' 2mapping/tmp.multi | \
		cut -d "p" -f 1 | \
		sort -u | \
		wc -l \
			> 2mapping/tmp.multi-loci.count
	
	# list of uce to remove
	awk '{print $1}' 2mapping/tmp.multi | \
		cut -d "p" -f 1 | \
		sort -u \
			> '2mapping/'$PROBE'_to_remove_2.list'

	# store summary
		if [ -f 2mapping/summary.txt ]; then

			printf "Writing mapping summary for $PROBE-$GEN_2\n"

		else

			printf "Writing mapping summary for $PROBE-$GEN_2\n"
		
			printf "probe-genome\tmapped-loci\ttotal-loci\tNumber-multi-mapped-loci\n" > 2mapping/summary_2.txt
		fi

	NAME=$(printf "$PROBE-$GEN")

	MAPPEDLOCI=$(cat 2mapping/tmp.uce.count)
	
	TOTALLOCI=$(cat '0data/probes/'$PROBE'.fasta' | grep ">" | cut -d "|" -f 1 | cut -d _ -f1 | sort -u | grep "." -c)
	
	MULTIMAPPING=$(cat 2mapping/tmp.multi-loci.count)

	printf "$NAME\t$MAPPEDLOCI\t$TOTALLOCI\t$MULTIMAPPING\n" >> 2mapping/summary_2.txt

	rm 2mapping/tmp.*

fi

printf "$PROBE-$GEN conversion and check complete\n\n"