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