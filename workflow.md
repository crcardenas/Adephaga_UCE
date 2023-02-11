# UCE-Characterization pipeline
Cody Raul Cardenas 02.2023

Annotated genomes from Darwin Tree of Life project
* GCA_944738965.1_Nebria_brevicollis
* GCA_943142095.1_Ophonus_ardosiacus
* GCA_911728475.2_Pterostichus_madidus

GFF files are available at: https://projects.ensembl.org/darwin-tree-of-life/

This following work flow assumes you have all files in the same be careful with the shell script, ensure you call conda properly. (friend mentioned something about mamba working more efficiently than conda)

**This script worked for both Nebria and Pterostichus, new tasks:**

1) print out example of compatable header for fasta file and gff file

2) take intput (e.g., probemap.sh -f $0 -g $1 -p $2 -o $3) where -f is fasta/fna, -g is the gff3 file, -p $2 probe file, and -o is the extension to save files as (e.g., nebriBrevi1).

3) get basic statistics about the total genes, total exons, the number of probes mapped to each respectively

4) get basic statistics about each UCE loci rather than each probe. (I can probably use bedtools to "merge" UCE matches)

---

Data needs formated as a first steps. The pipeline will need a genome.txt file for steps in bedtools. The fasta files and gff files do not have consistent headers so we need to format the fasta to match the GFF. For now I am going to stick with using shell, the resulting file will likely go to R or python for further statistics and UCE manipulation (concatenation at later steps).

naming convention uses the first four letters of genus and species (e.g., genuSpec#) and the genome version

---

## make new headers for FNA files
ensure your checking the scaffold/chromosome naming convention and adjust for each scaffold/chromosome awk replacement. May need to skip to the next step (make a bedtools genome file) in order to get naming convention
```
awk -F" " '{ if ($1~">CAJVRY") {gsub(">",""); print ">"$1} else {print $0}}' GCA_911728475.2_icPteMadi1.2_genomic.fna > GCA_911728475.2_icPteMadi1.2_genomic_new.fna

awk -F" " '{ if ($1~">OU") {print ">"$NF} else {print $0}}' GCA_911728475.2_icPteMadi1.2_genomic_new.fna > pterMadi2_newheader.fna
rm ./*_new.fna
```
## make a bedtools genome file
to get natsort use apt-get install python3-natsort, requires python3
This gets used later in bedtools
```
tr -s " " \\t < Pterostichus_madidus-GCA_911728475.1-2021_12-genes.gff | awk -F"\t" '{ if ($1~"##sequence-region") print $2 "\t" $4}' | natsort > pterMadi2.genomefile
```

## sort the gff file
here we need to identify the genomic features we are interested in our gff files to extract them from the base gff file. We are making two assumptions about the genome and genes:
1) intergenic sequences do not have any promotors close to genes
2) there is no alternative splicing or overlapping genes (ex overlaping genes: https://doi.org/10.1186/s12864-021-08039-6)

my conda environment contains bedtools, bwa, samtools
```
conda activate sam-bam-bedtools
bedtools sort -i Pterostichus_madidus-GCA_911728475.1-2021_12-genes.gff -g pterMadi2.genomefile > pterMadi2.sorted.gff
```
## extract gene features
We need to transform data from GFF (base 1) counting to bed counting (base 0), and print the gene feature ID (column 1, start -1, and end) using awk.
```
awk '$3 == "gene" {print $0}' pterMadi2.sorted.gff > pterMadi2.sorted.genes.gff
```

## excise intergenc sequences using bedtools complement
this function gets the complement ranges of gene features (e.g., whats not a gene is intergenic).
```
bedtools complement -i pterMadi2.sorted.genes.gff -g pterMadi2.genomefile > pterMadi2.sorted.intergenic.bed
```

## extract exons from gff file 
 need two files here, one with your typical bed files, and another with just exon in gff format. The gff format will allow for exon identifiers to be included in those gff files
```
awk '$3 == "exon" {print $1 "\t" $4-1 "\t" $5}' pterMadi2.sorted.gff > pterMadi2.sorted.exons.bed # this will be so we can map intergenic regions
awk '$3 == "exon" {print $0}' pterMadi2.sorted.gff > pterMadi2.sorted.exons.gff # this ensures that we have a gff file to give intersect file exon names for later workflow
```

## join and sort 
```
bedtools sort -i <(cat pterMadi2.sorted.intergenic.bed pterMadi2.sorted.exons.bed) -g pterMadi2.genomefile > pterMadi2.intergenic-and-exons.sorted.bed

```

## identify introns
get complements to identify regions that are introns using bedtools
```
bedtools complement -i pterMadi2.intergenic-and-exons.sorted.bed -g pterMadi2.genomefile > pterMadi2.sorted.introns.bed
```

---
# Identify UCEs that map to gene features
Now that we have our data formated and genomic characters seperated we can begin to map UCEs

## Map UCE probes to genome
Find where the UCE probes match in the genome using bwa. 
```
bwa index pterMadi2_newheader.fna
```

## alignt probes to genome using bwa-mem
```
bwa mem pterMadi2_newheader.fna Adephaga_2.9Kv1_UCE-Probes.fasta > pterMadi2-probes.sam
```

## convert sam to bam
samtools is already running in my environment
```
samtools view -h -b -S pterMadi2-probes.sam > pterMadi2-probes.mem.bam
```

## sort bam file
```
samtools sort pterMadi2-probes.mem.bam -o pterMadi2-probes.mem.sorted.bam
```
## extract sequenced that mapped against the genome
currently the bam file still contains unmapped information
```
samtools view -b -F 4 pterMadi2-probes.mem.sorted.bam > pterMadi2-probes.mem.sorted.mapped.bam
```
## create a bed file 
```
bedtools bamtobed -i pterMadi2-probes.mem.sorted.mapped.bam > pterMadi2-probes.mem.sorted.mapped.bed
```

## get intersect of adephaga probes and genome features
use the -names option in the same order as the -b files this will give us a generic annotation to use in later analyses
```
bedtools intersect -a pterMadi2-probes.mem.sorted.mapped.bed \
    -b pterMadi2.sorted.introns.bed \
    pterMadi2.sorted.exons.gff  \
    -names intron exon \
    -wb > Adephaga2.9-pterMadi2.introns-exons.intersect

bedtools intersect -a pterMadi2-probes.mem.sorted.mapped.bed \
    -b pterMadi2.sorted.intergenic.bed \
    pterMadi2.sorted.genes.gff \
    -names intergenic gene \
    -wb > Adephaga2.9-pterMadi2.intergenic-genentic.intersect
```

## format the intersect file
Columns are inconsistent otherwise, using a GFF file adds additional (but necessary) info. Here we get only what informatoin is necessary (not worried about strands or scores for now)
For every line that has ensembl print only those columns, all other lines should only have those columns


---

## intresect file format

```
scaffold    qstart    qend    query   type    seqname    seqstart   seqend  attribute
1	39310	39334	uce-127282_p11	60	+	intron	1	39281	39334
1	39334	39419	uce-127282_p11	60	+	exon	1	ensembl	exon	39335	39463	.	+	.	Parent=transcript:ENSMPTT00005003008;constitutive=0;exon_id=ENSMPTE00005013128;rank=8;version=1
1	39334	39419	uce-127282_p11	60	+	exon	1	ensembl	exon	39335	39463	.	+	.	Parent=transcript:ENSMPTT00005003000;constitutive=0;exon_id=ENSMPTE00005013128;rank=8;version=1
1	39339	39459	uce-127282_p12	60	+	exon	1	ensembl	exon	39335	39463	.	+	.	Parent=transcript:ENSMPTT00005003008;constitutive=0;exon_id=ENSMPTE00005013128;rank=8;version=1
1	39339	39459	uce-127282_p12	60	+	exon	1	ensembl	exon	39335	39463	.	+	.	Parent=transcript:ENSMPTT00005003000;constitutive=0;exon_id=ENSMPTE00005013128;rank=8;version=1
```

## general statistics using grep & wc in bash

_Pterostichus madidus_ and Adephaga 2.9k probes

total probes: 38948

total matches: 22072

exons: 11560 / 22072

intergenic: 7120 / 22072

introns: 3392 / 22072
