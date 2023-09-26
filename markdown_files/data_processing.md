# Data Processing

The `*.gff` and `*.fna` files from ensemble or genbank do not immediately cooperate with this pipeline and require data processing. This should be done in your 0data directory with those files.

Steps to process data:
1. clean up fna headers
2. sort `*.gff` files
3. create genome file
4. unrwap probe fasta files that need it

## new headers for FNA/FA/FASTA files

### `SEGMENTATION FAULT?`
```
index.sh: line 29: 72366 Segmentation fault      (core dumped)
```

BECAUSE YOU LOST ALL YOUR SEQUENCE DATA WITH YOU FASTA FILES WHEN EDITING THE HEADERS! CHECK THIS STEP AGAIN DINGUS!

check the scaffold/chromosome naming conventions in your `*.fna` and `*.gff3` files and adjust for replacement.

```
awk 'BEGIN { FS=" " } { if (/^>/) { print $1 } else { print $0 }}' genome.fa > gen1.fna
```


## bedtools genome file

```
tr -s " " \\t < Pterostichus_madidus-GCA_911728475.2-2022_03-genes.gff | awk -F"\t" '{ if ($1~"##sequence-region") print $2 "\t" $4}' | natsort > pterMadi2.genomefile
```

Alternatively `samtools faix` can be used.

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

# clean up tmp files
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