---
title: UCE characterization and concatenation
author: CR Cardenas
date: 2023.07
---

# TASKS

1) Run orthofinder with Spades data & probes
2) Run orthofinder with Spades data & new probeset
3) Incorporate distance within gene rather than generally distance. 
3) test&cleanup Rscript output
4) generate vendiagram of exons,introns,genic,intergenic loci as well (in R)
5) create json file with environment for download


# UCE characerization and concatenation for phylogenomic analysis of Adephaga

The goal is to map probes to genomes in order to realize the total overalp of all probes used n the Adephaga data for __concatenation/merging multiple UCEs__ (as in [Van Dam et al 2021](https://academic.oup.com/sysbio/article/70/2/307/5880562)) on the same gene and __UCE characterization__. In general, folks will only have one probeset to use. But because I am integrating [anchored hybrid enrichment data](https://resjournals.onlinelibrary.wiley.com/doi/full/10.1111/syen.12508), [Adephaga UCE probes](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5260)[Adephaga UCE probes](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5260), and the original [Coleoptera UCE probes](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12754) I will need to create a merged dataset (name TBD: cat/joined/frankenprobes).

Two important assumptions being made about the genome and genes being used:
1. intergenic sequences close to genes are not being considered as promotors
2. there is no alternative splicing or overlapping genes (ex: [overlaping genes](https://doi.org/10.1186/s12864-021-08039-6))


## Goals of this workflow:
1. identify what probes are genetic or intergenic 
2. identify the overlap/intersect of probesets used between datasets
3. create a new probeset that integrates all probes for use in [Phyluce](https://phyluce.readthedocs.io/en/latest/)
4. create a list of probes that should be concatenated in a partition file for phylogenetic analysis
    1. create a script to integrate it based on an existing partition file (e.g., output of Phyluce)


## Data

This pipeline will mainly use the [*Pterostichus madidus*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9315234/) genome for testing and finalizing the workflow, but others will be used. Genomes with annotations (general feature files, GFF) can be found at the [Darwin Tree of Life ensemble](https://projects.ensembl.org/darwin-tree-of-life/). 

Ensure you have a GFF3 format. For details on GFF format, see [here](http://gmod.org/wiki/GFF3).

### __*!!! genomes used are likely to change by final analysis, likely just the ones with GFF files Neb.brevicolis, Oph.ardosiacus, and Pte.madidus*__

__Available genomes and their names names__

|code|name|ID|
|---|---|---|
|agrpla1|*Agrilus planipennis*|GCA_000699045.1|
|\* amphizoa|*Amphizoa insolens*|DNA3784|
|anogla1|*Anoplophora glabripennis*|GCA_000390285.1|
|\* bemHap1|*Bembidion haplogonum*|DNA2544|
|\* chlSer1|*Chlaenius sericeus*|DNA4821|
|denpon1|*Dendroctonus ponderosae*|GCA_000355655.1|
|lepdec|*Leptinotarsa decemlineata*|GCA_000500325.1|
|lioTuu1|*Lionepha*|DNA3782|
|menmol1|*Mengenilla moldrzyki*|GCA_000281935.1|
|! nebBrevi1|*Nebria brevicolis*|GCA_944738965.1|
|! nebSali1|*Nebria salina*|GCA_944039245.1|
|\* omoHam1|*Omoglymmius hamatus*|DNA3783|
|onttau1|*Onthophagus taurus*|GCA_000648695.1|
|! ophArdo1|*Ophonus ardosiacus*|GCA_943142095.1|
|! pteMad2|*Pterostichus madidus*|GCA_911728475.2|
|\* pterMel1|*Pterostichus melenarius*|DNA3787|
|\* traGib1|*Trachypachus insolens*|DNA3786|
|tricas1|*Tribolium castaneum*|GCA_000002335.2|


(* indicates Adephaga reference genome; ! indicates a genome feature file)

__genomes:__
* GCA_002278615.1_Pchal_1.0_genomic.fna
* GCA_018344505.1_ASM1834450v1_genomic.fna
* GCA_002278615.1_Pchal_1.0_genomic.fna
* GCA_022063505.1_RBINS_Cgran_v1_genomic.fna
* GCA_911728475.2_icPteMadi1.2_genomic.fna
* GCA_913698125.1_Blapp1_genomic.fna
* GCA_913698245.1_Btrans1_genomic.fna
* GCA_933228885.1_icLeiSpin1.1_genomic.fna
* GCA_943142095.1_icOphArdo1.1_genomic.fna
* GCA_944039245.1_icNebSali1.1_genomic.fna
* GCA_944738965.1_icNebBrev1.1_genomic.fna
* GCA_947425015.1_icPteNige1.1_genomic.fna
* GCA_947534325.1_icAgoFuli1.1_genomic.fna

__probes:__
* vasilikopoulos_etal_ahe_probes.fasta
* Adephaga_2.9Kv1_UCE-Probes.fasta
* Coleoptera-UCE-1.1K-v1.fasta

### genomes are removed from the github repository due to their size

# Software and packages used

### __*!!! create a json file with environment for download and easy duplication of the environent*__

* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
* [BWA](https://bio-bwa.sourceforge.net/)
* [natsort](https://anaconda.org/anaconda/natsort)
* [sam tools](http://www.htslib.org/)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [R](https://www.r-project.org/)
* [OrthoFinder](https://github.com/davidemms/OrthoFinder)

Recomended installation procedure:
```
conda create --name characterization
conda activate characterization
conda install -c bioconda bedtools blast bwa samtools seqkit 
conda install -c anaconda natsort
```

This pipeline should be able to run on a personal laptop (windows with linux subsystem and linux, uncertain about mac) with sufficient storage, memory and CPU available. Alternatively, you can .... [add json information]

# Data Processing

The `*.gff` and `*.fna` files from ensemble or genbank do not immediately cooperate with this pipeline and require data processing. This should be done in your 0data directory with those files.

Steps to process data:
1. clean up fna headers
2. sort `*.gff` files
3. create genome file
4. unrwap probe fasta files that need it

## new headers for FNA files

check the scaffold/chromosome naming conventions in your `*.fna` files and adjust for awk replacement outlined here:

```
awk -F" " '{ if ($1~">CAJVRY") {gsub(">",""); print ">"$1} else {print $0}}' GCA_911728475.2_icPteMadi1.2_genomic.fna > tmp.fna


awk -F" " '{ if ($1~">OU") {print ">"$NF} else {print $0}}' tmp.fna > pteMadi2.fna


rm ./tmp.fna
```

## bedtools genome file


```
tr -s " " \\t < Pterostichus_madidus-GCA_911728475.2-2022_03-genes.gff | awk -F"\t" '{ if ($1~"##sequence-region") print $2 "\t" $4}' | natsort > pterMadi2.genomefile
```

## sort GFF file
```
bedtools sort -i Pterostichus_madidus-GCA_911728475.2-2022_03-genes.gff -g pterMadi2.genomefile > pterMadi2.sorted.gff
```

## process probe files

some probe files are wrapped, this pipeline assumes you have an unwrapped fasta file (see coleoptera UCE)
```
seqkit seq -w 0 Coleoptera-UCE-1.1K-v1.fasta > Coleoptera-UCE-1.1k.fasta
```

to simplify the pipeline and naming conventions I am shortening the probe files as well

## directories structure:

```
.
├── 0data
│   ├── genomes
│   │   ├── pteMadi2.fna
│   │   ├── pteMadi2.genomefile
│   │   ├── pteMadi2.sorted.gff
│   │   ├── Pterostichus_madidus-GCA_911728475.2-2022_03-genes.gff
│   │   └── Pterostichus_madidus-GCA_911728475.2-softmasked.fa
│   └── probes
│       ├── Adephaga_2.9Kv1_UCE-Probes.fasta
│       ├── Adephaga.fasta
│       ├── Coleoptera.fasta
│       ├── Coleoptera-UCE-1.1k.fasta
│       ├── Coleoptera-UCE-1.1K-v1.fasta
│       ├── Vasil2021_AHE_probes.fasta
│       └── Vasil.fasta
└── workflow.md
```


!!!! once pipeline is finalized get genomes with GFF files into the genome directory !!!!

## Extract gene features

The `*.gff` files do not have all genome features of interest so they need to be identified, this can be done with simple `awk` and `bedtools` commands.

This script takes the users `*.gff` file, uses the input name to seperate genic, intergenic, exon, and intron features from the sorted `*.gff` file. Ensure your naming convention is consistent because it looks for that file in your `./0data/genomes` directory not in your currnet working diectory. (e.g., if your input `*.gff` file is named `pteMadi2.sorted.gff` your sorted.gff file pteMadi2 will be added too all output files based on the first field).

To run the following script create a directory 1genefeatures in the same directory as 0data:
```
bash genefeatures.sh pteMadi1.sorted.gff
```
this script can also be run with a for loop with a list if you have multiple gff files:
```
for i in $(cat mylist.txt); do 
	bash genefeatures.sh $i
done
```

### `genefeatures.sh`

```
#!/bin/bash
# take user input
GFFINPUT=$1
GFFNAME=$(echo $GFFINPUT | cut -d "." -f 1)

#activate environment
source /local/anaconda3/bin/activate
conda activate characterization

# extract gene features using awk into gff and bed format
printf "extracting gene features of ${GFFINPUT}\n"
awk '$3 == "gene" {print $0}' '0data/genomes/'${GFFNAME}'.sorted.gff' \
	> '1genefeatures/'${GFFNAME}'.sorted.genes.gff'
awk '$3 == "gene" {print $1 "\t" $4-1 "\t" $5 "\t""gene"}' '1genefeatures/'${GFFNAME}'.sorted.genes.gff' \
	> '1genefeatures/'${GFFNAME}'.sorted.genes.bed'

# extract intergenic features
printf "extracting intergenic features of ${GFFINPUT}\n"
bedtools complement -i '1genefeatures/'${GFFNAME}'.sorted.genes.gff' \
	-g '0data/genomes/'${GFFNAME}'.genomefile' \
	> '1genefeatures/'${GFFNAME}'.sorted.intergenic.tmp' # this line gives coordinates for all intergenic regions
awk '{print $0 "\t" "intergenic"}' '1genefeatures/'${GFFNAME}'.sorted.intergenic.tmp' \
	> '1genefeatures/'${GFFNAME}'.sorted.intergenic.bed' # create a bedfile of just intergenic regions
awk '{print $1 "\t" "bedtools_complement" "\t" "intergenic" "\t" $2+1 "\t" $3 "\t" "." "\t" "." "\t" "." "\t" "description=intergenic-region-IDd-with-bedtools"}' '1genefeatures/'${GFFNAME}'.sorted.intergenic.bed' \
	> '1genefeatures/'${GFFNAME}'.sorted.intergenic.gff' # create gff file of intergenic regions

# extract exon features
printf "extracting exon features of ${GFFINPUT}\n"
awk '$3 == "exon" {print $1 "\t" $4-1 "\t" $5 "\t" "exon"}' '0data/genomes/'${GFFNAME}'.sorted.gff' \
	> '1genefeatures/'${GFFNAME}'.sorted.exons.bed' # generate a simple exon bed file, watch the count! GFF is 1-count and BED is 0-count
awk '$3 == "exon" {print $0}' '0data/genomes/'${GFFNAME}'.sorted.gff' \
	> '1genefeatures/'${GFFNAME}'.sorted.exons.gff' # this ensures that we have a gff file to give intersect file exon names for a later step

# join and sort intergenic and exonic regions to extract intron
printf "extracting intron features of ${GFFINPUT}\n"
bedtools sort -i <(cat '1genefeatures/'${GFFNAME}'.sorted.intergenic.bed' '1genefeatures/'${GFFNAME}'.sorted.exons.bed') \
	-g '0data/genomes/'${GFFNAME}'.genomefile' \
	> '1genefeatures/'${GFFNAME}'.sorted.intergenic-and-exons.bed.tmp' # this sorts the input of the concatenated features, allows introns to be identified with bedtools
bedtools complement -i '1genefeatures/'${GFFNAME}'.sorted.intergenic-and-exons.bed.tmp' \
	-g '0data/genomes/'${GFFNAME}'.genomefile' \
	> '1genefeatures/'${GFFNAME}'.sorted.introns.bed' # create complement bed file
awk '{print $1 "\t" "bedtools_complement" "\t" "intron" "\t" $2+1 "\t" $3 "\t" "." "\t" "." "\t" "." "\t" "description=intron-region-IDd-with-bedtools"}' '1genefeatures/'${GFFNAME}'.sorted.introns.bed' \
	> '1genefeatures/'${GFFNAME}'.sorted.introns.gff' # create gff file of intronic features

# create new gff files
printf "creating new GFF file for ${GFFNAME}\n"
cat '1genefeatures/'${GFFNAME}'.sorted.introns.gff' \
	'1genefeatures/'${GFFNAME}'.sorted.exons.gff' \
	'1genefeatures/'${GFFNAME}'.sorted.intergenic.gff' \
	'1genefeatures/'${GFFNAME}'.sorted.genes.gff' \
	> '1genefeatures/'${GFFNAME}'_simple.gff.tmp'
bedtools sort -i '1genefeatures/'${GFFNAME}'_simple.gff.tmp' \
	-g '0data/genomes/'${GFFNAME}'.genomefile' \
	> '1genefeatures/'${GFFNAME}'_simple.sorted.gff'

cat '0data/genomes/'${GFFNAME}'.sorted.gff' \
	'1genefeatures/'${GFFNAME}'.sorted.introns.gff' \
	> '1genefeatures/'${GFFNAME}'_detailed.gff.tmp'
bedtools sort -i '1genefeatures/'${GFFNAME}'_detailed.gff.tmp' \
	-g '0data/genomes/'${GFFNAME}'.genomefile' \
	> '1genefeatures/'${GFFNAME}'_detailed.sorted.gff'

# clean up tmp filescardenas
rm 1genefeatures/*.tmp
printf "output files in 1genefeatures:\n\t ${GFFNAME}_simple.sorted.gff is a new gff file with just intergenic, genic, exon, and intron information\n\t ${GFFNAME}_detailed.sorted.gff is a new gff file with intron inforation as well as original CDS, mrna, and 3' + 5' UTR info\n\n"
```
* Edit this shell script to properlly call your environment. On my server using `$ source genefeatures.sh` or `$ bash -i genefeatures.sh` will still run this script. But it may throw a conda error, I'm not sure how to fix this. Keep this in mind if you run into errors with further scripts.

### Current directory structure

```
.
├── 0data
│   ├── genomes
│   │   ├── pteMadi2.fna
│   │   ├── pteMadi2.genomefile
|	|	└── ...
│   └── probes
│       ├── Adephaga_2.9Kv1_UCE-Probes.fasta
│       ├── Adephaga.fasta
|		└── ...
├── 1genefeatures
│   ├── pteMadi2_detailed.sorted.gff
│   ├── pteMadi2_simple.sorted.gff
|	└── ...
├── genefeatures.sh
└── workflow.md
```

# Map probes to gene features and intersect

### __*!!! integrate blast into the `index.sh` script*__

Depending on the detail needed for downstream analysis or description you can use the simpe or detailed gff files to map your probes too. For simplicity, I will be using the `*_simple.sorted.gff` file. There will be a few shell scripts that talk to each other that produce what probes intersect with gene features in the genome as well as what probes dont intersect. 

To run this script first __make two directories `2mapping` and `3intersect`__ in the same directory as `0data` and `1genefeatures`

You can use a single genome and your probes of interest:

```
intersect.sh genome.fna probes1.fasta probes2.fasta probes3.fasta
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

This script first indexes the genome and maps using BWA and then does the same using using blastn.

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

The indexed `*.sam` is converted to a bed file and their is a check for probes with multiple positions in the genome; potential paralogs. The output then returns a final `*.bed` file, summary of mapping, and list of those probes/loci which should be excluded based on the sam file. The script creates a new probeset fasta file if there are probes (and the entire loci). and re-maps them 


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

	printf "output $PROBE-$GEN.mem.sorted.mapped.bed in 2mapping/ directory\n"

fi

if [ -f 2mapping/*.bam ]; then

rm 2mapping/*.bam

printf "Checking for potential loci/probbes map to different regions in $PROBE-$GEN.sam\n"

else

printf "Checking for potential loci/probbes map to different regions in $PROBE-$GEN.sam\n"

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
awk '{print $1}' 2mapping/tmp.multi | cut -d "_" -f 1 | sort -u | wc -l > 2mapping/tmp.multi-uce.count

# list of uce to remove
awk '{print $1}' 2mapping/tmp.multi | cut -d "_" -f 1 | sort -u > '2mapping/'$PROBE'_to_remove.list'

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

	printf "Loci removed fom $PROBE\t" && grep ">" '0data/probbes/'$PROBE'.fasta' -c

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
	
		printf "Checking if loci/probbes map to different regions in $PROBE-$GEN.sam\n"
	
	else
	
		printf "Checking if loci/probbes map to different regions in $PROBE-$GEN.sam\n"
	
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
		cut -d "_" -f 1 | \
		sort -u | \
		wc -l \
			> 2mapping/tmp.multi-loci.count
	
	# list of uce to remove
	awk '{print $1}' 2mapping/tmp.multi | \
		cut -d "_" -f 1 | \
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
	bash intersect.sh $i \
	probes1-genome_2.mem.sorted.mapped.bed \
	probes2-genome_2.mem.sorted.mapped.bed \
	probes3-genome.mem.sorted.mapped.bed;
done
```

The example here is set up for 3 probe files as input, but you could alternatively set it up to work with more or less (see commented out lines in the script). The order of your `*.bed` files is important here if using more than one. Use the probeset of most interest first, as any additional probes will be compared to this first probeset. In the example above, Adephaga is first. 

Lastly, the mapped loci are checked for probes that do and do not intersect in the genome using `bedtools`. The output returns an `*.intersect` file for all probes compared to the first input probeset, and then for each additional combination that do not intersect (e.g. a `nointersect` file for: probe1-probe2-probe3, probe2-probe1-probe3, probe3-probe1-probe2).

```
source /local/anaconda3/bin/activate
conda activate characterization

# turn off N probes (e.g., PROBE1:PROBEN) depending on probesets used for characterization

# user bedfile input
BED1=$2
BED2=$3
BED3=$4
# create variable for naming probes and genomes
GEN=$(echo $1 | cut -d . -f 1)
PROBE1=$(echo $2 | cut -d - -f 1)
PROBE2=$(echo $3 | cut -d - -f 1)
PROBE3=$(echo $4 | cut -d - -f 1)


printf "getting intersect of all probes with probe1\n"
# all probes intersect
bedtools intersect \
	-a '2mapping/'$BED1 \
	-b '2mapping/'$BED2 \
	'2mapping/'$BED3 \
	-names $PROBE2 $PROBE3 \
	-wa -wb > '3intersect/'$GEN'-PROBES.intersect';
# turn off N -b inputs depending on the total number of probes you use

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

Ensure intersect files have a consistent format. Gene feature information should be retained for features that have it (name, ID, parent transcript, etc.). Here is an example script 

```
awk '{ if ($9 =="ensembl") {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$11"\t"$12"\t"$16} else {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$10"\t"}}' Adephaga-pteMadi2.introns-exons.intersect > Adephaga2.9-pterMadi2.introns-exons.out.intersect

awk '{ if ($9 =="ensembl") {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$11"\t"$12"\t"$16} else {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$10"\t"}}' Adephaga2.9-pterMadi2.intergenic-genentic.intersect > Adephaga2.9-pterMadi2.intergenic-genentic.out.intersect
```

Generate a list that identifies the proportion of probes mapped and the total genes:

```
awk '$3 == "gene" {print $1 "\t" $4-1 "\t" $5-1 "\t" $3 "\t" $9}' pterMadi2.sorted.gff > pterMadi2.sorted.all_genes.gff
```


## Current Directory Structure
The directory structure should now look like this 

```
.
├── 0data
│   ├── genomes
│   │   ├── pteMadi2.fna
│   │   ├── pteMadi2.fna.amb
│   │   ├── pteMadi2.fna.ann
│   │   └── ...
│   └── probes
│       ├── Adephaga_2.9Kv1_UCE-Probes.fasta
│       ├── Adephaga_2.fasta
│       ├── Adephaga.fasta
│       └── ...
├── 1genefeatures
│   ├── pteMadi2_detailed.sorted.gff
│   ├── pteMadi2_simple.sorted.gff
│   ├── pteMadi2.sorted.exons.bed
│   └── ...
├── 2mapping
│   ├── Adephaga-pteMadi2_2.mem.sorted.mapped.bed
│   ├── Adephaga-pteMadi2_2.sam
│   ├── Adephaga-pteMadi2.mem.sorted.mapped.bed
│   └── ...
├── 3intersect
│   ├── pteMadi2-Coleoptera.nointersect
│   ├── pteMadi2-Coleoptera-nointersect-summary.txt
│   ├── pteMadi2-nointersect-summary.txt
│   └── ...
├── genefeatures.sh
├── index.sh
├── intersect.sh
├── run.sh
├── sam2bed.sh
├── summary.txt
└── workflow.md
```


...

# Integrate probesets

If you find it necessary to join multiple probesets follow the rest of this pipeline. Some steps of this pipeline will provide a summary of the probes that may be useful.

## Subset overlap between probefiles

### `blast.sh`

order is important
if you have an probefile that isn't in the expected phyluce format (e.g., `uce100_p10 | ...`) then you need to modify the cut in the last awk statements or do this step manually. This file takes input [....]

```
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
	cut -d _ -f 1 | \
	sort -u \
	> '4probe-reduction/'$PROBEFILE2'.to-remove'

PROBE2NAME2=$(echo $2 | cut -d _ -f 1)

fi

```

Next add to the `*.to-remove` list based on intersect files. In this example I am prioritizing Adephaga probes.

```
awk '$7 == "Coleoptera" {print $11}' 3intersect/pteMadi2-PROBES.intersect | cut -d _ -f 1 | sort -u >> 4probe-reduction/Coleoptera_2.to-remove
awk '{print $0}' 4probe-reduction/Coleoptera_2.to-remove | sort -u > 4probe-reduction/Coleoptera_2.to-remove.list 
```

In this example, the Vasil probes are formated differently and for downstream purposes I am considering exons over the genes (note the change in cut)

```
awk '$7 == "Vasil" {print $11}' 3intersect/pteMadi2-PROBES.intersect | cut -d _ -f 1,2,3 | sort -u >> 4probe-reduction/Vasil.to-remove
awk '{print $0}' 4probe-reduction/Vasil.to-remove | sort -u > 4probe-reduction/Vasil.to-remove.list
```

next we remove the duplicates using the generated lists

```
#cat 0data/probes/Adephaga_2.fasta 0data/probes/Coleoptera_2.fasta 0data/probes/Vasil.fasta > 4probe-reduction/tmp.fasta
#cat 4probe-reduction/*.to-remove.list > 4probe-reduction/tmp.list

awk '{ if ((NR>1) && ($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' 0data/probes/Coleoptera_2.fasta | grep -vf 4probe-reduction/Coleoptera_2.to-remove.list - | tr "\t" "\n" > 4probe-reduction/Coleoptera.subset.fasta

awk '{ if ((NR>1) && ($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' 0data/probes/Vasil.fasta | grep -vf 4probe-reduction/Vasil.to-remove.list - | tr "\t" "\n" > 4probe-reduction/Vasil.subset.fasta
```

## rename probe files 

In order to use these probes for the phyluce pipeline, the probe ID's need renamed for non-adephaga probesets (Coleoptera and Vasil). First the original name is appended to the probe ID at the end of the header.

```
# ensure you have a new directory 5join-probes

awk 'BEGIN{RS=">";-F" |"} (NR>1) {print ">" $1 " " $2",original-probe-ID:"$1 "\n" $3}' 4probe-reduction/Coleoptera.subset.fasta > 5join-probes/Coleoptera.subset.tmp

awk 'BEGIN{RS=">";-F" |"} (NR>1) {print ">" $1 " " $2",original-probe-ID:"$1 "\n" $3}' 4probe-reduction/Vasil.subset.fasta > 5join-probes/Vasil.subset.tmp

grep ">" 5join-probes/Coleoptera.subset.tmp | cut -d"|" -f1 | tr -d ">" > 5join-probes/Coleoptera-rename.list

grep ">" 5join-probes/Vasil.subset.tmp | cut -d"|" -f1 | tr -d ">" > 5join-probes/Vasil-rename.list
```

Next a delimited file is created, with the old probe ID, new target loci, new probe number, and final renamed loci. Recall that the Vasil probes, TC######-PA_# is the "gene" and TC######-PA_#\_# exon, and the rest, TC######-PA_#\_#\_#\_#_# probe ID's. This may vary depending on the probes used.

!!!! If using microsoft visual studio code then you should check the first line of the resulting files: __the first line may have strange input inserted into it.__ !!!!

```
for i in $(cat 5join-probes/Coleoptera-rename.list); do
EXON=$(echo $i | cut -d _ -f 1);
echo -e $i "\t" $EXON;
done > 5join-probes/Coleoptera-rename.tmp.list

for i in $(cat 5join-probes/Vasil-rename.list); do
LOCI=$(echo $i | cut -d _ -f 1,2,3);
echo -e $i "\t" $LOCI;
done > 5join-probes/Vasil-rename.tmp.list
```

Next create a new ID for our probes based on their "loci" so it is compatable with phyluce. This uses the new tmp.list file we created previously to alter the fasta file.

```
awk 'BEGIN {count=1000000;probe[0]=""};
# if the probe is not loaded into the array; that probe gets +1
{if (!($2 in probe)) {probe[$2]=count++} ;
# print a new column with your new loci name
print $0, "\t", "uce-"count}' 5join-probes/Coleoptera-rename.tmp.list > 5join-probes/Coleoptera-rename.tmp.list2

# do the same, but for our newly created column
# compare the value from the first field stored in save
awk '$3 != save { counter = 1; save = $3 }
# if the values differ resets to one
{ print $0, "\t", "_p"counter++ }' 5join-probes/Coleoptera-rename.tmp.list2 > 5join-probes/Coleoptera-rename.tmp.list3

# and kiss
awk '{print ">"$1, "\t", ">"$3 $4}' 5join-probes/Coleoptera-rename.tmp.list3 > 5join-probes/Coleoptera.renamed-probe.list

#again for Vasil but at 2 million so there is no overlap
awk 'BEGIN {count=2000000;probe[0]=""};
{if(!($2 in probe)) {probe[$2]=count++};
print $0, "\t", "uce-"count}' 5join-probes/Vasil-rename.tmp.list > 5join-probes/Vasil-rename.tmp.list2
awk '$3 != save { counter = 1; save = $3 }
{ print $0, "\t", "_p" counter++ }' 5join-probes/Vasil-rename.tmp.list2 > 5join-probes/Vasil-rename.tmp.list3
awk '{print ">" $1, "\t", ">" $3 $4}' 5join-probes/Vasil-rename.tmp.list3 > 5join-probes/Vasil.renamed-probe.list

# clean up tmp files
rm 5join-probes/*.tmp.list*
```
!!!! if using microsoft visual studio code then you should check the first line of the resulting files: __the first line may have strange input inserted into it.__ !!!!

Lastly we concatenate the probe files

```
cat 0data 0data/probes/Adephaga_2.fasta 5join-probes/Coleoptera.subset.tmp 5join-probes/Vasil.subset.tmp > 5join-probes/joined_probes.tmp.fasta
cat 5join-probes/Coleoptera.renamed-probe.list 5join-probes/Vasil.renamed-probe.list > 5join-probes/joined.renamed-probes.list

awk 'BEGIN {OFS="\t"} NR==FNR {a[$1]=$2; next} /^>/ {if (a[$1]) {$1=a[$1]}}1' 5join-probes/joined.renamed-probes.list 5join-probes/joined_probes.tmp.fasta | tr "\t" " " > 5join-probes/joined.renamed.tmp.probes
```

following this pipeline the resulting probeset contains 4387 loci comprised of 91675 probes. So the new probefile will follow the naming convention of other UCE probesets joined_probes4.3k.fasta

```
mv 5join-probes/joined.renamed.tmp.probes > 5join-probes/joined_probes4.3k.fasta
```

Here are the resulting changes for the probes

`$ seqkit stats *.fasta`

|file|format|type|num_seqs|sum_len|min_len|avg_len|max_len|
|---|---|---|---|---|---|---|---|
|0data/probes/Adephaga_2.9Kv1_UCE-Probes.fasta|FASTA|DNA|38,948|4,673,760|120|120|120|
|0data/probes/Adephaga_2.fasta|FASTA|DNA|38,409|4,609,080|120|120|120|
|0data/probes/Adephaga.fasta|FASTA|DNA|38,948|4,673,760|120|120|120|
|0data/probes/Coleoptera_2.fasta|FASTA|DNA|13,559|1,627,080|120|120|120|
|0data/probes/Coleoptera.fasta|FASTA|DNA|13,674|1,640,880|120|120|120|
|0data/probes/Coleoptera-UCE-1.1k.fasta|FASTA|DNA|13,674|1,640,880|120|120|120|
|0data/probes/Vasil2021_AHE_probes.fasta|FASTA|DNA|49,786|5,974,320|120|120|120|
|0data/probes/Vasil.fasta|FASTA|DNA|49,786|5,974,320|120|120|120|
|5join-probes/Coleoptera.subset.tmp|FASTA|DNA|7,303|876,360|120|120|120|
|5join-probes/Vasil.subset.tmp|FASTA|DNA|45,963|5,515,560|120|120|120|
|5join-probes/joined_probes4.3k.fasta|FASTA|DNA|91,675|11,001,000|120|120|120|


### Current directory structure

```
.
├── 0data
│   ├── genomes
│   │   ├── pteMadi2.fna
│   │   ├── pteMadi2.fna.amb
│   │   ├── pteMadi2.fna.ann
│   │   └── ...
│   └── probes
│       ├── Adephaga_2.9Kv1_UCE-Probes.fasta
│       ├── Adephaga_2.fasta
│       ├── Adephaga_2.fasta.ndb
│       └── ...
├── 1genefeatures
│   ├── pteMadi2_detailed.sorted.gff
│   ├── pteMadi2_simple.sorted.gff
│   ├── pteMadi2.sorted.exons.bed
│   └── ...
├── 2mapping
│   ├── Adephaga-pteMadi2_2.mem.sorted.mapped.bed
│   ├── Adephaga-pteMadi2_2.sam
│   ├── Adephaga-pteMadi2.mem.sorted.mapped.bed
│   └── ...
├── 3intersect
│   ├── pteMadi2-Coleoptera.nointersect
│   ├── pteMadi2-Coleoptera-nointersect-summary.txt
│   ├── pteMadi2-nointersect-summary.txt
│   └── ...
├── 4probe-reduction
│   ├── Coleoptera_2-Adephaga_2.blastn
│   ├── Coleoptera_2.to-remove
│   ├── Coleoptera_2.to-remove.list
│   └── ...
├── 5join-probes
│   ├── Coleoptera.renamed-probe.list
│   ├── Coleoptera-rename.list
│   ├── Coleoptera.subset.tmp
│   └── ...
├── blast.sh
├── genefeatures.sh
├── index.sh
├── intersect.sh
├── run.sh
├── sam2bed.sh
└── workflow.md
```

# Check orthology of data
## Jeremy's idea

use UCE & Vasil probes along side the assembled spades loci to check for orthogroups, we want to keep UCe loci (uce & vasil exons) that have unique uces (e.g., only orhtogroups contianing UCE1) not paralogs (UCE10 & 11 in the same orthogroup). With the 300 coleoptera-adephaga overlap.

## Overview
the goal here is to ensure that we have selected orthologous probes. Adephaga probes have 300 Coleoptera loci, so there should be approximately a 300 shared uces. I was only able to find ~ 200, this may be due to the fidelity of the Coleoptera probekit OR it may be due to a poor pipeline. 

28 all spades contigs from those that previously had a high number of recovered loci to `orthofinder/contigs` directory and included the trimmed probes

need to include sequence info in the header, otherwise its nonsense!

` for f in *.fasta; do sed -i "s/^>/>${f%.*}_/" "$f"; done `

` orthofinder -f contigs2/ -t 10 -a 2 -og -d `

```
##~#~#~#~#~#
#26.04.2023#
##~#~#~#~#~#
conversation with grey
check just Adephaga probes... new trimming: clipkit, relaxed gblocks trimming, bossert found a way to tweak overlap

check out untrimmed alignments in just Adephaga2.9k kit
pull loci using probes, treat dataset outside of phyluce, relaxed gblocks, gblocks, try clipkit, tweak mafft parameters think about reducing overlap in AHE probes... 
AHE probes are pulling loci 

Van Dam et al (2021) used taxa specific probes, but here we are using a broad number of taxa to ensure that overlap is consistent across taxa, not just one genome.
```


# Concatenate Loci for phylogenetic analysis

Two choices, concatenate loci that are found on the same gene OR concatenate loci based on their distance from one another. The assumption for by gene is that loci on the gene should have a similar evolutionary rate. However, genes can have varying rates of evolution within them. Creating a partition by length, assuming linkage between loci, it is possible to join loci that need joined.

** notes **

aught to remap the newprobeset *OR* rename the probes before mapping. It makes sense to remap here. 

mkdir 6concatenate-loci

1) map and make bedfile 
2) find clusters using `bedtools cluster -d` and check a range from 0:1000 and examine the change in the number of "clusters"
3) write script to concatenate loci/uce's by the cluster statistic chosen

```
cat 2mapping/Adephaga-pteMadi2_2.mem.sorted.mapped.bed 2mapping/Coleoptera-pteMadi2_2.mem.sorted.mapped.bed 2mapping/Vasil-pteMadi2.mem.sorted.mapped.bed | bedtools sort | bedtools cluster -d 1000  > tmp.cluster1000

cat 2mapping/Adephaga-pteMadi2_2.mem.sorted.mapped.bed 2mapping/Coleoptera-pteMadi2_2.mem.sorted.mapped.bed 2mapping/Vasil-pteMadi2.mem.sorted.mapped.bed | bedtools sort | bedtools cluster -d 500  > tmp.cluster500
```

# Concatenate loci by gene


### `group-by-gene.sh`

This code takes in user input and maps the new probes against the genome. It produces a list of loci that group together on the same gene. For ease, it is best to simplify the created probefile to be succinct (joined_probes4.3k.fasta > joined.fasta)

```
bash group-by.sh pteMadi2.fna joined.fasta
```

```
#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

GEN=$(echo $1 | cut -d . -f 1)
PROBE=$(echo $2 | cut -d . -f 1)

# index genome and map probes
if [ -f '6group-by/'###'sam' ]; then

	printf "$PROBE already mapped to -$GEN \n\n";

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
		
	printf "finished bwa indexing and mapping\noutput $PROBE-$GEN.sam in 6group-by directory \n\n";

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
    -wa -wb > '6group-by/'${PROBE}-${GEN}'.gene.intersect'

#here we group the GFF feature info and sumarize by UCE loci
# -g group by column or ranges (1-4; 1,2,3,4)
# -c summarize by this column
printf "group $PROBE mapped to $GEN \n"
bedtools groupby -i '6group-by/'${PROBE}-${GEN}'.gene.intersect' \
	 -g 9 -c 13 -o collapse > '6group-by/'${PROBE}'-'${GEN}'.grpby_gene'

# create output files

cat '6group-by/'${PROBE}'-'${GEN}'.grpby_gene' | tr -d " " > 6group-by/tmp
while read TMP; do
	GENE=$(echo $TMP | cut -d";" -f 1 | cut -d: -f 2)
    UCELIST=$(echo $TMP | awk '{print $2}')
    echo -e "$GENE \t $UCELIST"
done < 6group-by/tmp | tr -d " " > 6group-by/tmp2

GENOME=$(echo $GEN | cut -d"_" -f1,2)
awk 'BEGIN{FS=OFS="\t"} 																		# set intput and output as tab delimited
	{split($2, probes, ","); 																	# split the second feild of each line into an array called probes, and as a comma delimited file
	for (i in probes) {split(probes[i], loci, "_"); 											# loop through each element in the proe array, and split it according to the underscore
	loci_array[loci[1]] = 1} delete probes; 													# while in the for loop save the first part of the loci as associative array called loci_array
	for (j in loci_array) {loci_list = loci_list ? loci_list "," j : j}; 						# loop through each key in the loci array and append a string called "loci list". The string is built up using the ? : in the loop condition: if loci_list is not empty (if it does not contain a previous loci) the current loci is appended with a comma; otherwise the current loci is added to the string.
	delete loci_array; 																			# delete the loci array to free up memory for the next line
	print "'$GENOME'", $1, loci_list; loci_list=""}' 6group-by/tmp2 > '6group-by/'${PROBE}'-'${GEN}'.short.grpby_gene'		# print the first field ($1) followed by the loci_list string built earlier and then reset loci_list to start with the next line.

awk 'BEGIN{FS=OFS="\t"} {split($2, probes, ","); for (i in probes) {split(probes[i], loci, "_"); loci_array[loci[1]] = 1} delete probes; for (j in loci_array) {loci_list = loci_list ? loci_list "," j : j}; delete loci_array; print "'$GENOME'", $1, loci_list; loci_list=""}'

printf "finished grouping $PROBE by $GEN genes \n"
printf "\t Detailed grouped by file is ${PROBE}-${GEN}.grpby_gene\n"
printf "\t Succinct grouped by gene file is ${PROBE}-${GEN}.short.grpby_gene\n"
```
this gives us a filethat looks like 

```
ENSMPTG0000500001   uce-1,uce-2,uce-3,uce-7
ENSMPTG0000700105   uce-4,uce-5,uce-6
ENSMPTG0001400123   uce-9,uce-10
...
```

Using a nexus partition file like this:
```
#nexus
begin sets;
charset uce-1.nexus = 1-100;
charset uce-2.nexus = 101-200;
charset uce-3.nexus = 201-300;
charset uce-4.nexus = 301-400;
charset uce-5.nexus = 401-500;
charset uce-6.nexus = 500-1500;
charset uce-7.nexus = 1501-1700;
charset uce-8.nexus = 1700-10000;
charset uce-9.nexus = 10001-10100;
charset uce-10.nexus = 10100-10500;
done;
```

We can run a script that generates a new partition of our sequence data:

```
#nexus
begin sets;
charset ENSMPTG0000500001 = 1-100 101-200 201-300 1501-1700;
charset ENSMPTG0000700105 = 301-400 401-500 500-1500;
charset uce-8 = 1700-10000;
charset ENSMPTG0001400123 = 10001-10100 10100-10500;
done;
```

First make a directory to work from 7modify-partition, and copy your phyluce partition file into it. Because phyluce uses the file name, and not just the loci/uce name, it needs the `*.nexus` removed from the file. Additionally, we only want to rename UCEs if there is more than one found on a gene. Notably, depending on the completeness of the matrix that gets put in, then many taxa might be missing one or more of this loci/probe 

```
cat 7modify-partition/inttrim_30p_parti.nex | sed 's/.nexus//g' > 7modify-partition/tmp.nexus
grep "," 6group-by/joined-pteMadi2.short.grpby_gene > 7modify-partition/tmp.grpby_gene
```

### `modify-partitions.awk`

to run this script use the `*.short.grpby_gene`

```
awk -f modify-partitions.awk 7modify-partition/tmp.grpby_gene 7modify-partition/tmp.nexus > 7modify-partition/new.nexus
```

```
#!/usr/bin/awk -f


# Read the gene mapping file into an array
FNR == NR {
    split($2, uces, ",");
    for (i in uces) {
        gene_map[uces[i]] = $1;
    }
    next;
}

# Process the partition file
/^charset/ {
    charset_name = $2;
    charset_pos = substr($0, index($0, "=") + 2);
    if (charset_name in gene_map) {
        gene = gene_map[charset_name];
        if (gene in new_charset) {
            new_charset[gene] = new_charset[gene] " " charset_pos;
        } else {
            new_charset[gene] = charset_pos;
        }
    } else {
        new_charset[charset_name] = charset_pos;
    }
}

END {
    print "#nexus";
    print "begin sets;";
    for (gene in new_charset) {
        gsub(";", "", new_charset[gene]);
        printf "charset %s = %s;\n", gene, new_charset[gene];
    }
    print "done;";
}
```

# Concatenate by distance

### `group-by-cluster.sh`

1) need to find intergenic loci that are within N base pairs from one another
2) need to find genic loci that are within N base pairs from one another
3) need to exclude any loci that are genic and intergenic

In this workflow we need to map by both genes and intergenic regions so we do not cluster both genic and intergenic regions
We have already created a genic mapping, but the intergenic regions needs the same treatment.
The first input is your genome, the second is the name of the probe fasta file and the third is the window for clustering (e.g., within 500 bp of the UCEs)
```
bash group-by-cluster.sh pteMadi2.fna joined.fasta 500
```

```
#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

GEN=$(echo $1 | cut -d . -f 1)
PROBE=$(echo $2 | cut -d . -f 1)
BP=$3

# get intergenic probes (we already have the intersect of genic probes)
echo -e "Find intersect of $PROBE and $GEN intergenic features \n"
bedtools intersect -a '1genefeatures/'${GEN}'.sorted.intergenic.gff' \
    -b '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bed' \
    -names ${PROBE} \
    -wb > '6group-by/'${PROBE}-${GEN}'.intergenic.intersect'

# get bed files of intergenic and genic loci
awk -F"\t" 'BEGIN{OFS="\t"} {print $10,$11,$12,$13,$14,$15}' '6group-by/'${PROBE}'-'${GEN}'.gene.intersect' > 6group-by/tmp.gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} {print $10,$11,$12,$13,$14,$15}' '6group-by/'${PROBE}'-'${GEN}'.intergenic.intersect' > 6group-by/tmp.intergenic.bed

# find genic and intergenic regions that intersect
bedtools intersect -a 6group-by/tmp.intergenic.bed \
	-b 6group-by/tmp.gene.bed | 
	awk '{print $4}' | \
	cut -d _ -f 1 | \
	sort -u \
	> 6group-by/genic-intergenic_loci.list

# get cluster from bedfile
bedtools cluster -i '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bed' \
	-d $BP > 6group-by/probes.cluster

# exclude intergenic-genic loci from cluster list
grep -vf 6group-by/genic-intergenic_loci.list 6group-by/probes.cluster > 6group-by/probes.filtered.cluster

# get bedtools cluster file grouped by the cluster
bedtools groupby -g 7 -c 4 -o collapse \
	-i 6group-by/probes.filtered.cluster \
	> 6group-by/tmp.probes.grpby_cluster

awk 'BEGIN{FS=OFS="\t"} {split($2, probes, ","); 
	for (i in probes) {split(probes[i], loci, "_"); 
	loci_array[loci[1]] = 1} delete probes; 
	for (j in loci_array)  {loci_list = loci_list ? loci_list "," j : j}; 
	delete loci_array; 
	print "cluster_"$1, loci_list; loci_list=""}' 6group-by/tmp.probes.grpby_cluster > '6group-by/'${PROBE}'-'${GEN}'.grpby_cluster'

```

Just like before, the only loci of interest are those that cluster together within 500bp

```
grep "," 6group-by/joined-pteMadi2.short.grpby_cluster > 7modify-partition/tmp.grpby_cluster
```

### `modify-partitions.awk`

since the loci of interest have already been removed, the awk script can be run to create a partition file where loci within 500 bp of each other are partitioned together

```
awk -f modify-partitions.awk 7modify-partition/tmp.grpby_cluster 7modify-partition/tmp.nexus > 7modify-partition/cluster.nexus
```

```
#!/usr/bin/awk -f


# Read the gene mapping file into an array
FNR == NR {
    split($2, uces, ",");
    for (i in uces) {
        gene_map[uces[i]] = $1;
    }
    next;
}

# Process the partition file
/^charset/ {
    charset_name = $2;
    charset_pos = substr($0, index($0, "=") + 2);
    if (charset_name in gene_map) {
        gene = gene_map[charset_name];
        if (gene in new_charset) {
            new_charset[gene] = new_charset[gene] " " charset_pos;
        } else {
            new_charset[gene] = charset_pos;
        }
    } else {
        new_charset[charset_name] = charset_pos;
    }
}

END {
    print "#nexus";
    print "begin sets;";
    for (gene in new_charset) {
        gsub(";", "", new_charset[gene]);
        printf "charset %s = %s;\n", gene, new_charset[gene];
    }
    print "done;";
}
```


### System and Environment Information

```
Linux pyrgus 4.4.0-119-generic #143-Ubuntu SMP Mon Apr 2 16:08:24 UTC 2018 x86_64 x86_64 x86_64 GNU/Linux

     active environment : characterization
    active env location : /home/cody/.conda/envs/characterization
            shell level : 2
       user config file : /home/cody/.condarc
 populated config files : /home/cody/.condarc
          conda version : 4.6.14
    conda-build version : 3.17.6
         python version : 3.7.1.final.0
       base environment : /local/anaconda3  (read only)
           channel URLs : https://conda.anaconda.org/conda-forge/linux-64
                          https://conda.anaconda.org/conda-forge/noarch
                          https://conda.anaconda.org/bioconda/linux-64
                          https://conda.anaconda.org/bioconda/noarch
                          https://repo.anaconda.com/pkgs/main/linux-64
                          https://repo.anaconda.com/pkgs/main/noarch
                          https://repo.anaconda.com/pkgs/free/linux-64
                          https://repo.anaconda.com/pkgs/free/noarch
                          https://repo.anaconda.com/pkgs/r/linux-64
                          https://repo.anaconda.com/pkgs/r/noarch
          package cache : /local/anaconda3/pkgs
                          /home/cody/.conda/pkgs
       envs directories : /home/cody/.conda/envs
                          /local/anaconda3/envs
               platform : linux-64
             user-agent : conda/4.6.14 requests/2.21.0 CPython/3.7.1 Linux/4.4.0-119-generic ubuntu/16.04.4 glibc/2.23
                UID:GID : 1028:1040
             netrc file : None
           offline mode : False

# packages in environment at /home/cody/.conda/envs/characterization:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
bedtools                  2.31.0               hf5e1c6e_2    bioconda
blast                     2.14.0          pl5321h6f7f691_2    bioconda
bwa                       0.7.17              he4a0461_11    bioconda
bzip2                     1.0.8                h7f98852_4    conda-forge
c-ares                    1.18.1               h7f98852_0    conda-forge
ca-certificates           2023.01.10           h06a4308_0    anaconda
certifi                   2022.12.7                pypi_0    pypi
curl                      7.87.0               h5eee18b_0  
entrez-direct             16.2                 he881be0_1    bioconda
gettext                   0.21.1               h27087fc_0    conda-forge
htslib                    1.6                  h42e7767_0    anaconda
keyutils                  1.6.1                h166bdaf_0    conda-forge
krb5                      1.19.3               h3790be6_0    conda-forge
ld_impl_linux-64          2.38                 h1181459_1    anaconda
libcurl                   7.87.0               h91b91d3_0  
libdeflate                1.14                 h166bdaf_0    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libev                     4.33                 h516909a_1    conda-forge
libffi                    3.4.2                h6a678d5_6    anaconda
libgcc-ng                 12.2.0              h65d4601_19    conda-forge
libgomp                   12.2.0              h65d4601_19    conda-forge
libidn2                   2.3.4                h166bdaf_0    conda-forge
libnghttp2                1.47.0               hdcd2b5c_1    conda-forge
libssh2                   1.10.0               haa6b8db_3    conda-forge
libstdcxx-ng              12.2.0              h46fd767_19    conda-forge
libunistring              0.9.10               h7f98852_0    conda-forge
libuuid                   1.41.5               h5eee18b_0    anaconda
libzlib                   1.2.13               h166bdaf_4    conda-forge
natsort                   7.1.1              pyhd3eb1b0_0    anaconda
ncbi-vdb                  2.11.0               h1b792b2_1    bioconda
ncurses                   6.3                  h27087fc_1    conda-forge
openssl                   1.1.1s               h7f8727e_0    anaconda
pcre                      8.45                 h9c3ff4c_0    conda-forge
perl                      5.26.2            h36c2ea0_1008    conda-forge
perl-archive-tar          2.32                    pl526_0    bioconda
perl-carp                 1.38                    pl526_3    bioconda
perl-common-sense         3.74                    pl526_2    bioconda
perl-compress-raw-bzip2   2.087           pl526he1b5a44_0    bioconda
perl-compress-raw-zlib    2.087           pl526hc9558a2_0    bioconda
perl-exporter             5.72                    pl526_1    bioconda
perl-extutils-makemaker   7.36                    pl526_1    bioconda
perl-io-compress          2.087           pl526he1b5a44_0    bioconda
perl-io-zlib              1.10                    pl526_2    bioconda
perl-json                 4.02                    pl526_0    bioconda
perl-json-xs              2.34                    pl526_1    bioconda
perl-list-moreutils       0.15                    pl526_1    bioconda
perl-pathtools            3.75            pl526h14c3975_1    bioconda
perl-scalar-list-utils    1.52            pl526h516909a_0    bioconda
pip                       22.3.1                   pypi_0    pypi
python                    3.10.8               h7a1cb2a_1    anaconda
readline                  8.2                  h5eee18b_0    anaconda
samtools                  1.6                  hcd7b337_9    bioconda
seqkit                    2.5.0                h9ee0642_0    bioconda
setuptools                65.6.3                   pypi_0    pypi
sqlite                    3.40.1               h5082296_0    anaconda
tk                        8.6.12               h1ccaba5_0    anaconda
tzdata                    2022a                hda174b7_0    anaconda
wget                      1.20.3               ha56f1ee_1    conda-forge
wheel                     0.37.1             pyhd3eb1b0_0    anaconda
xz                        5.2.6                h166bdaf_0    conda-forge
zlib                      1.2.13               h166bdaf_4    conda-forge

```