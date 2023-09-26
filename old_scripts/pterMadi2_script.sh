#!/bin/bash

##### 
# new workflow. 
# Need to extract just pterostichus from the probes
# clunky way to clean them up for R analysis:
awk -F"| " '{ if ($1 ~ ">") {gsub(",", "|"); print $0} else {print $0}}' Adephaga_2.9Kv1_UCE-Probes.fasta | awk -F"|" '{if ($1 ~ ">") {print $1 "|" $6} else {print $0}}' |  awk -F "|" '{if ($1 ~ ">") {gsub("probes-source:", ""); print $1 "|" $2} else {print $0}}' > Adephaga2.9k_probes_short.fasta
    # that will let us do some quick stats in r to calculate percentage of spp mapped to Pterostichus genome
awk -F "|" '{ if ($2 == "pterMel") {print $1 "|" $2 NR%1}}'  Adephaga2.9k_probes_short.fasta 

# subset by source 
# this requires that you have only 
awk -v i=pterMel1 ' BEGIN{j="^>.*" i "[0-9]*"} $0 ~ j {k=1; print; next} /^>/ {k=0} k {print} ' Adephaga_2.9Kv1_UCE-Probes.fasta > pterMel1_short.fasta


# genomic source names:
#   ensure youre probe file has the same number of delimiters!!!!
# awk -F "," '!/^$/ {print $5}' Adephaga_2.9Kv1_UCE-Probes.fasta | cut -d":" -f 2| sort -u
# (*) indicates original Adephaga probes, GCA_* come from https://doi.org/10.1111/2041-210X.12754
# agrpla1   Agrilus planipennis GCA_000699045.1
# * amphizoa  Amphizoa insolens DNA3784
# anogla1   Anoplophora glabripennis    GCA_000390285.1
# * bemHap1   Bembidion haplogonum  DNA2544
# * chlSer1   Chlaenius sericeus    DNA4821
# denpon1   Dendroctonus ponderosae GCA_000355655.1
# lepdec1   Leptinotarsa decemlineata   GCA_000500325.1
# lioTuu1   Lionepha    DNA3782
# menmol1   Mengenilla moldrzyki    GCA_000281935.1
# * omoHam1   Omoglymmius hamatus   DNA3783
# onttau1   Onthophagus taurus  GCA_000648695.1
# * pterMel1  Pterostichus melenarius   DNA3787
# * traGib1   Trachypachus insolens DNA3786
# tricas1   Tribolium castaneum   GCA_000002335.2

# need to get the percentage of probes that mapped to the Pterostichus genome and what species
# they were designed from

## EXPECTED IN SAME DIRECTORY ##
# probe file
# genomicfasta/scaffoldfasta
# annotated_gff file

# subset UCE loci from pterostichus spp


# make new header
awk -F" " '{ if ($1~">CAJVRY") {gsub(">",""); print ">"$1} else {print $0}}' GCA_911728475.2_icPteMadi1.2_genomic.fna > GCA_911728475.2_icPteMadi1.2_genomic_new.fna

awk -F" " '{ if ($1~">OU") {print ">"$NF} else {print $0}}' GCA_911728475.2_icPteMadi1.2_genomic_new.fna > pterMadi2_newheader.fna
rm ./*_new.fna

#get genome file !!! need to ensure you have natsort installed on your system/environment
tr -s " " \\t < Pterostichus_madidus-GCA_911728475.2-2022_03-genes.gff | awk -F"\t" '{ if ($1~"##sequence-region") print $2 "\t" $4}' | natsort > pterMadi2.genomefile


# sort GFF file
#source /home/nepenthes/miniconda3/etc/profile.d/conda.sh
#conda activate sam-bam-bedtools
source /home.nepenthes/mambaforge/etc/profile.d/conda.sh 
mamba init
conda activate sam-bam-bedtools
bedtools sort -i Pterostichus_madidus-GCA_911728475.2-2022_03-genes.gff -g pterMadi2.genomefile > pterMadi2.sorted.gff

# generate intergenic regions
#awk '$3 == "gene" {print $1 "\t" $4-1 "\t" $5}' pterMadi2.sorted.gff > pterMadi2.sorted.genes.bed
awk '$3 == "gene" {print $0}' pterMadi2.sorted.gff > pterMadi2.sorted.genes.gff

# excise intergenic regions from gene features
bedtools complement -i pterMadi2.sorted.genes.gff -g pterMadi2.genomefile > pterMadi2.sorted.intergenic.bed

# extract exons
# need two files here, one with your typical bed files, and another with just exon in gff format. The gff format will allow for exon identifiers to be included in those gff files
awk '$3 == "exon" {print $1 "\t" $4-1 "\t" $5-1}' pterMadi2.sorted.gff > pterMadi2.sorted.exons.bed # this will be so we can map intergenic regions
awk '$3 == "exon" {print $0}' pterMadi2.sorted.gff > pterMadi2.sorted.exons.gff # this ensures that we have a gff file to give intersect file exon names for later workflow

# join and sort
bedtools sort -i <(cat pterMadi2.sorted.intergenic.bed pterMadi2.sorted.exons.bed) -g pterMadi2.genomefile > pterMadi2.intergenic-and-exons.sorted.bed

# excise intronic features
bedtools complement -i pterMadi2.intergenic-and-exons.sorted.bed -g pterMadi2.genomefile > pterMadi2.sorted.introns.bed

# now to mapping probes
# index genome
bwa index pterMadi2_newheader.fna

# align probes to indexed genome
bwa mem pterMadi2_newheader.fna Adephaga_pterMele1.fasta > probes.sam

# convert sam to bam
samtools view -h -b -S probes.sam > probes.mem.bam

# sort bam file
samtools sort probes.mem.bam -o probes.mem.sorted.bam

# extract sequences that map against the genome
samtools view -b -F 4 probes.mem.sorted.bam > probes.mem.sorted.mapped.bam

# create bed file
bedtools bamtobed -i probes.mem.sorted.mapped.bam > probes.mem.sorted.mapped.bed

#! Problems with R transformations of GFF attribute column, need to make two seperate files
# mutate them, and then rejoin for our final file!
# get intersect of probes and gene features(introns or exons)
bedtools intersect -a probes.mem.sorted.mapped.bed \
    -b pterMadi2.sorted.introns.bed \
    pterMadi2.sorted.exons.gff  \
    -names intron exon \
    -wb > Adephaga2.9-pterMadi2.introns-exons.intersect

# get intersect of probes as intergenic or genenic
bedtools intersect -a probes.mem.sorted.mapped.bed \
    -b pterMadi2.sorted.intergenic.bed \
    pterMadi2.sorted.genes.gff \
    -names intergenic gene \
    -wb > Adephaga2.9-pterMadi2.intergenic-genentic.intersect

# need to format the intersect file
# if not columns are consistent, using a GFF file adds additional info. Here we get only what informatoin is necessary (not worried about strands or scores for now)
# For every line that has ensembl print only those columns, all other lines should only have those columns
awk '{ if ($9 =="ensembl") {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$11"\t"$12"\t"$16} else {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$10"\t"}}' Adephaga2.9-pterMadi2.introns-exons.intersect > Adephaga2.9-pterMadi2.introns-exons.out.intersect
awk '{ if ($9 =="ensembl") {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$11"\t"$12"\t"$16} else {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$10"\t"}}' Adephaga2.9-pterMadi2.intergenic-genentic.intersect > Adephaga2.9-pterMadi2.intergenic-genentic.out.intersect
