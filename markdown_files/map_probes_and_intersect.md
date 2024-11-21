# Map probes to gene features and intersect

Depending on the detail needed for downstream analysis or description you can use the simpe or detailed gff files to map your probes too. For simplicity, I will be using the `*_simple.sorted.gff` file. There will be a few shell scripts that talk to each other that produce what probes intersect with gene features in the genome as well as what probes dont intersect. 

To run this script first __make two directories `2mapping` and `3intersect`__ in the same directory as `0data` and `1genefeatures`

You can use a single genome and your probes of interest:

```
bash run.sh genome.fna probes1.fasta probes2.fasta probes3.fasta
```

Or you can use a list of the genomes you are interested in.Importantly, you need to modify `run.sh` to include more than two probes!
```
for i in $(cat genomes\.list); do 
	bash run.sh $i \
	probes1.fasta \
	probes2.fasta \
	probes3.fasta;
done
```


### `run.sh`
This script runs the pipeline, calling the `index.sh` and then `sam2bed.sh` scripts. These files should be in the same directory as the `run.sh`. Ensure you input the genome fna file 

```
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
```

### `index.sh`

This script first indexes the genome and maps using BWA.

```
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
	printf "mapping $PROBE to $GEN \n";
	bwa mem '0data/genomes/'$GEN'.fna' '0data/probes/'$PROBE'.fasta' \
		> '2mapping/bwa/'$PROBE'-'$GEN'.sam';

else
	printf "indexing $GEN \n";
	bwa index '0data/genomes/'$GEN'.fna';
	
	printf "\n mapping $PROBE to $GEN \n";
	bwa mem '0data/genomes/'$GEN'.fna' '0data/probes/'$PROBE'.fasta' \
		> '2mapping/'$PROBE'-'$GEN'.sam';
	
	printf " finished bwa indexing, output $PROBE-$GEN.sam in 2mapping directory \n\n";
fi

mkdir blast


```

### `sam2bed.sh`

The indexed `*.sam` is converted to a bed file and their is a check for probes with multiple positions in the genome; potential paralogs. The output then returns a final `*.bed` file, summary of mapping, and list of those probes/loci which should be excluded based on the sam file. The script creates a new probeset fasta file if there are probes (and the entire loci) to remove and re-maps them.


```
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

	printf "output $PROBE-$GEN.mem.sorted.mapped.bed in 2mapping/ directory\n\n"

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
	cut -d "p" -f 1 | \
	sort -u | \
	grep "." -c \
		> '2mapping/'tmp.uce.count

printf "Total number of $PROBE loci mapped against $GEN\t" && cat 2mapping/tmp.uce.count

# number of multiple mappings
# count multi mapping
cat '2mapping/'$PROBE'-'$GEN'.sam' | 
	grep -v "@" | \
	grep "SA:" -c \
		> 2mapping/tmp.multi-uce.count

printf "Total $PROBE probes with multiple mappings in $GEN\t" && cat 2mapping/tmp.multi-uce.count

# store  names
cat '2mapping/'$PROBE'-'$GEN'.sam' | grep -v "@" | grep "SA:" > 2mapping/tmp.multi

# count UCE/Loci names
awk '{print $1}' 2mapping/tmp.multi | cut -d "p" -f 1 | sort -u | wc -l > 2mapping/tmp.multi-uce.count

# list of uce to remove
awk '{print $1}' 2mapping/tmp.multi | cut -d "p" -f 1 | sort -u > '2mapping/'$GEN'-'$PROBE'_to_remove.list'

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
	grep -vf '2mapping/'$GEN'-'$PROBE'_to_remove.list' - | \
	tr "\t" "\n" \
		> '0data/probes/'$PROBE'_2.fasta';

# grep -A 1 -vf '2mapping/'$PROBE'_to_remove.list' '0data/probes/'$PROBE'.fasta' > '0data/probes/'$PROBE'_2.fasta'

	printf "Probes removed from $PROBE\t" && grep "_" '2mapping/'$GEN'-'$PROBE'_to_remove.list' -c

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
		cut -d "p" -f 1 | \
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

			printf "Writing mapping summary for $PROBE-$GEN""_2\n"

		else

			printf "Writing mapping summary for $PROBE-$GEN""_2\n"
		
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

exit
```

### `intersect.sh`

The `intersectn.sh` expects an `*.fna` genome file and the `*.bed` files that were generated. Ensure you call the `*_2.mem.sorted.mapped.bed` file if the probeset was altered.

e.g.,
```
bash intersect.sh genome.fna probes1-genome_2.mem.sorted.mapped.bed probes2-genome_2.mem.sorted.mapped.bed probes3-genome.mem.sorted.mapped.bed
```

Again, multiple genomes can be used if you want

```
for i in $(cat genomes\.list); do
	GEN=$(echo $i | cut -d "." -f 1)
	bash intersect.sh $i \
	2mapping/Adephaga-${GEN}_2.mem.sorted.mapped.bed \
	2mapping/Coleoptera-${GEN}_2.mem.sorted.mapped.bed \
	2mapping/Vasil-${GEN}.mem.sorted.mapped.bed;
done
```

The example here is set up for 3 probe files as input, but you could alternatively set it up to work with more or less (see commented out lines in the script). The order of your `*.bed` files is important here if using more than one. Use the probeset of most interest first, as any additional probes will be compared to this first probeset. In the example above, Adephaga is first. 

Lastly, the mapped loci are checked for probes that do and do not intersect in the genome using `bedtools`. The output returns an `*.intersect` file for all probes compared to the first input probeset, and then for each additional combination that do not intersect (e.g. a `nointersect` file for: probe1-probe2-probe3, probe2-probe1-probe3, probe3-probe1-probe2).

```
#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

# turn off N probes (e.g., PROBE1:PROBEN) depending on probesets used for characterization

# user bedfile input
BED1=$2
BED2=$3
BED3=$4
#BED#=$#

# create variable for naming probes and genomes
GEN=$(echo $1 | cut -d . -f 1)
PROBE1=$(echo $2 | cut -d . -f 1 | cut -d - -f 1)
PROBE2=$(echo $3 | cut -d . -f 1 | cut -d - -f 1)
PROBE3=$(echo $4 | cut -d . -f 1 | cut -d - -f 1)
#PROBE#=$(echo $# | cut -d . -f 1 | cut -d - -f 1)

printf "getting intersect of all probes with probe1\n"
# all probes intersect
bedtools intersect \
	-a  '1genefeatures/'$GEN'_simple.sorted.gff' \
	-b '2mapping/'$BED1 \
	'2mapping/'$BED2 \
	'2mapping/'$BED3 \
	-names $PROBE1 $PROBE2 $PROBE3 \
	-wa -wb > '3intersect/'$GEN'-PROBES.intersect';
# turn off the number of -b inputs depending on the total number of probes you use
#	e.g., using only two your names might look like: -names $PROBE1 $PROBE2 \
# summary
bedtools summary -i '3intersect/'$GEN'-PROBES.intersect' \
	-g '0data/genomes/'$GEN'.genomefile' \
	> '3intersect/'$GEN'-PROBES-intersect-summary.txt'

printf "getting which probes do not intersect with probe1\n"
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

printf "getting which probes do not intersect with probe2\n"
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

printf "getting which probes do not intersect with probe3\n"
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
```

## Format intersect files

You can reduce redunancy in the `genome-PROBES.intersect` file given that probes mapped to genes will be reported twice, as genic and either exon or introns. You can create two files one with intergenic and genic mappings and one with intergenic, introns, and exons mappings. However, these can be filtered simply using R. 

It may be necessary to change what is included and excluded based on whether the `genome_detailed.sorted.gff` or `genome_simple.sorted.gff` file was used. This example code is shown works with the `*_simple.sorted.gff` and will produce a an tab delimited file with all probefiles. To subset by a condition simply select based on the column: 

```
awk -F"\t" '($3 != "gene") {print $11"\t"$12"\t"$13"\t"$14"\t"$3"\t"$1"\t"$4"\t"$5"\t""probeset="$10";"$9}' file.intersect | awk -F";" '$1=$1' OFS="\t" > out.intersect
```

The intersect of the probe flies can now be analyized with the R script or the rest of the pipline can be followed
