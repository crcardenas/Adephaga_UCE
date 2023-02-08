#!/bin/bash

# tasks for this script
# 1 - print out example of compatable header for fasta file and gff file
# 2 - take intput (e.g., probemap.sh -f $0 -g $1 -p $2 -o $3) where -f is fasta/fna, -g is the gff3 file, -p $2 probe file, and -o is the extension to save files as (e.g., nebriBrevi1).
# 3 - make basic statistics about the total genes, total exons, the number of probes mapped to each respectively

# make new header
#awk -F" " '{ if ($1~">CAJVRY") {gsub(">",""); print ">"$1} else {print $0}}' GCA_911728475.2_icPteMadi1.2_genomic.fna > GCA_911728475.2_icPteMadi1.2_genomic_new.fna

#awk -F" " '{ if ($1~">OU") {print ">"$NF} else {print $0}}' GCA_911728475.2_icPteMadi1.2_genomic_new.fna > pterMadi2_newheader.fna
#rm ./*_new.fna

#get genome file
#tr -s " " \\t < Pterostichus_madidus-GCA_911728475.1-2021_12-genes.gff3 | awk -F"\t" '{ if ($1~"##sequence-region") print $2 "\t" $4}' | natsort > pterMadi2.genomefile

# previous tasks completed before creation of this script

# sort GFF file
source /home/nepenthes/miniconda3/etc/profile.d/conda.sh
conda activate sam-bam-bedtools
bedtools sort -i Pterostichus_madidus-GCA_911728475.1-2021_12-genes.gff -g pterMadi2.genomefile > pterMadi2.sorted.gff

# generate intergenic regions
#awk '$3 == "gene" {print $1 "\t" $4-1 "\t" $5}' pterMadi2.sorted.gff > pterMadi2.sorted.genes.bed
awk '$3 == "gene" {print $0}' pterMadi2.sorted.gff > pterMadi2.sorted.genes.gff

# excise intergenic regions from gene features
bedtools complement -i pterMadi2.sorted.genes.gff -g pterMadi2.genomefile > pterMadi2.sorted.intergenic.bed

# extract exons
# need two files here, one with your typical bed files, and another with just exon in gff format. The gff format will allow for exon identifiers to be included in those gff files
awk '$3 == "exon" {print $1 "\t" $4-1 "\t" $5}' pterMadi2.sorted.gff > pterMadi2.sorted.exons.bed # this will be so we can map intergenic regions
awk '$3 == "exon" {print $0}' pterMadi2.sorted.gff > pterMadi2.sorted.exons.gff # this ensures that we have a gff file to give intersect file exon names for later workflow

# join and sort
bedtools sort -i <(cat pterMadi2.sorted.intergenic.bed pterMadi2.sorted.exons.bed) -g pterMadi2.genomefile > pterMadi2.intergenic-and-exons.sorted.bed

# excise intronic features
bedtools complement -i pterMadi2.intergenic-and-exons.sorted.bed -g pterMadi2.genomefile > pterMadi2.sorted.introns.bed

# now to mapping probes
# index genome
bwa index pterMadi2_newheader.fna

# align probes to indexed genome
bwa mem pterMadi2_newheader.fna Adephaga_2.9Kv1_UCE-Probes.fasta > pterMadi2-probes.sam

# convert sam to bam
samtools view -h -b -S pterMadi2-probes.sam > pterMadi2-probes.mem.bam

# sort bam file
samtools sort pterMadi2-probes.mem.bam -o pterMadi2-probes.mem.sorted.bam

# extract sequences that map against the genome
samtools view -b -F 4 pterMadi2-probes.mem.sorted.bam > pterMadi2-probes.mem.sorted.mapped.bam

# create bed file
bedtools bamtobed -i pterMadi2-probes.mem.sorted.mapped.bam > pterMadi2-probes.mem.sorted.mapped.bed

# get intersect of probes and genomic features
bedtools intersect -a pterMadi2-probes.mem.sorted.mapped.bed \
    -b pterMadi2.sorted.intergenic.bed \
    pterMadi2.sorted.introns.bed \
    pterMadi2.sorted.exons.gff  \
    -names intergenic intron exon \
    -wb > Adephaga2.9-pterMadi2.intersect
