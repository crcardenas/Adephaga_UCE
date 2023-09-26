---
title: Concatenated probes for integrative phylogenomics
author: CR Cardenas
date: 2023 05 17
---

# Generate non-intersecting probes

Need to map multiple genomes to realize the overalp of all probesets. 

Probe files need to be unwrapped, a single line (use seqkit) for later processing it is useful to have youre fna & gff files have the same name for chromosomes/scaffolds as well as having genome files. genome files describe the length of each chromosome/scaffold and the name of that chromosome/scaffold, see [bedtools](https://bedtools.readthedocs.io/en/latest/). 

The AHE probes have the following format: 
TC007324-PA_1 is the gene TC007324-PA_1_4 is the exon, all other codes are probe identifiers

`vasilikopoulos_etal_ahe_probes.fasta` was generated using their probefile set.

**probes:**		
- vasilikopoulos_etal_ahe_probes.fasta
- Adephaga_2.9Kv1_UCE-Probes.fasta
- Coleoptera-UCE-1.1K-v1.fasta

**genomes:**
- GCA_002278615.1_Pchal_1.0_genomic.fna
- GCA_018344505.1_ASM1834450v1_genomic.fna
- GCA_022063505.1_RBINS_Cgran_v1_genomic.fna
- GCA_911728475.2_icPteMadi1.2_genomic.fna
- GCA_913698125.1_Blapp1_genomic.fna
- GCA_913698245.1_Btrans1_genomic.fna
- GCA_933228885.1_icLeiSpin1.1_genomic.fna
- GCA_943142095.1_icOphArdo1.1_genomic.fna
- GCA_944039245.1_icNebSali1.1_genomic.fna
- GCA_944738965.1_icNebBrev1.1_genomic.fna
- GCA_947425015.1_icPteNige1.1_genomic.fna
- GCA_947534325.1_icAgoFuli1.1_genomic.fna

code requires the following software: [samtools](http://www.htslib.org/), [bedtools](https://bedtools.readthedocs.io/en/latest/), and [bwa](https://bio-bwa.sourceforge.net/index.shtml). 
useful tools used but unnecssary unless you perform downstream analyses: [seqkit](https://bioinf.shenwei.me/seqkit/)

I have been using mambaforge as my package mamanger, so use at your own discression.

---

# bash scripts used

to run the intersect script for one use just the genome, and probe files you are trying to intersect 

```
intersect.sh genome.fa probes1.fasta probes2.fasta probes3.fasta
```

or for multiple genomes,and what I ran: use a list of multiple genomes:

```
for i in $(cat genomes\.list); do 
	bash intersect.sh $i Adephaga_2.9Kv1_UCE-Probes.fasta vasilikopoulos_etal_ahe_probes.fasta Coleoptera-UCE-1.1K-v1.fasta
done
```

### ` intersect.sh `

```
!/bin/bash
source  ~/mambaforge/etc/profile.d/conda.sh
conda activate sam-bam-bedtools
GEN=$1
PROBE1=$2
PROBE2=$3

bash control.sh $GEN $PROBE1;
bash control.sh $GEN $PROBE2;

bedtools intersect \
	-a $PROBE1-$GEN.mem.sorted.mapped.bed \
	-b $PROBE2-$GEN.mem.sorted.mapped.bed $PROBE3-$GEN.mem.sorted.mapped.bed\
	-names $PROBE2 $PROBE3 \
	-wa -wb > $GEN.probes.intersect;

bedtools intersect \
	-a $PROBE2-$GEN.mem.sorted.mapped.bed \
	-b $PROBE1-$GEN.mem.sorted.mapped.bed $PROBE3-$GEN.mem.sorted.mapped.bed\
	-names $PROBE1 $PROBE3 \
	-v > $GEN.$PROBE2.nointersect;

bedtools intersect \
	-a $PROBE3-$GEN.mem.sorted.mapped.bed \
	-b $PROBE2-$GEN.mem.sorted.mapped.bed $PROBE1-$GEN.mem.sorted.mapped.bed\
	-names $PROBE2 $PROBE1 \
	-v > $GEN.$PROBE3.nointersect;
```

### ` control.sh `
```
!/bin/bash
source  ~/mambaforge/etc/profile.d/conda.sh
conda activate sam-bam-bedtools
GEN=$1
PROBE=$2

if [ -f $GEN.pac ]; then  
    bwa mem $GEN $PROBE > $PROBE-$GEN.sam;
    else
    bash index.sh $GEN;
    bwa mem $GEN $PROBE > $PROBE-$GEN.sam;
fi
```


### ` index.sh `

```
#!/bin/bash
source  ~/mambaforge/etc/profile.d/conda.sh
conda activate sam-bam-bedtools
GEN=$1
bwa index $GEN
```


### ` convert.sh `

```
!/bin/bash
source  ~/mambaforge/etc/profile.d/conda.sh
conda activate sam-bam-bedtools
GEN=$1
PROBE=$2

# convert sam to bam
samtools view -h -b -S $PROBE-$GEN.sam > $PROBE-$GEN.mem.bam

# sort bam file
samtools sort $PROBE-$GEN.mem.bam -o $PROBE-$GEN.mem.sorted.bam

# extract sequences that map against the genome
samtools view -b -F 4 $PROBE-$GEN.mem.sorted.bam > $PROBE-$GEN.mem.sorted.mapped.bam

# create bed file
bedtools bamtobed -i $PROBE-$GEN.mem.sorted.mapped.bam > $PROBE-$GEN.mem.sorted.mapped.bed
```

The resulting intersect files are what probes target the same regions in the genome used

Additionally, you get what probes *dont* intersect with the same regions as the UCE data. With this data we can identify what probes we should concatenate into a new franken-probe set by concatenating all nointersect and filtering by unique probes.

Should consider that we want to include ALL probes from that targeted loci/gene not just partial matches. Duplicates will be removed by phyluce

first unwrap fasta files that need it (multi-line to single line fasta)

` seqkit seq -w 0 Coleoptera-UCE-1.1K-v1.fasta > Coleoptera-UCE-1.1K-v1.tmp.fasta `

now use the list of non-intersecting loci to generate probeset

```
for i in $(cat Coleoptera-UCE-1.1K-v1.fasta.nointersect); do 
awk '/'"$i"'/{print;getline;print}' Coleoptera-UCE-1.1K-v1.tmp.fasta; 
done > Coleoptera-UCE-1.1k_nointersect.fasta

for i in $(cat vasilikopoulos_etal_ahe_probes.fasta.nointersect); do 
awk '/'"$i"'/{print;getline;print}' vasilikopoulos_etal_ahe_probes.fasta; 
done > vasilikopoulos_etal_ahe_nointersect.fasta
```

join all probes together!
```
cat vasilikopoulos_etal_ahe_nointersect.fasta Adephaga_2.9Kv1_UCE-Probes.fasta Coleoptera-UCE-1.1k_nointersect.fasta > franken_probes.fasta
```

There *MAY* be issues integrating this custom probeset into phyluce but we wont know until we do it!

something strange in the formatting here... need to double check our fasta files careful how you use cat and your awk script... seomthing caused at least 2 lines that I could find to paste the header at the end of a fasta file...

---

# Generate unique identifiers

phyluce requires a particular format for each probe to properly use the pipeline
` uce-###_p## `

Need to identify the max ID number in Adephaga probes, then jump through the millions to adequately seperate probe numbers. We can append the real name to the  metadata after the pipe
e.g.: 
``` 
>uce-500_p1 |design:coleoptera-v1,designer:faircloth,probes-locus:uce-500,probes-probe:1,probes-source:lepdec1,probes-global-chromo:KI578650.1,probes-global-start:149735,probes-global-end:149855,probes-local-start:10,probes-local-end:130 
actgactg 
> uce...
```

turns into 

``` 
>uce-500000000_p1 |design:coleoptera-v1,designer:faircloth,probes-locus:uce-500,probes-probe:1,probes-source:lepdec1,probes-global-chromo:KI578650.1,probes-global-start:149735,probes-global-end:149855,probes-local-start:10,probes-local-end:130,original-probe:uce-500_p1
actgactg 
> uce...
```
This allows us to keep track of the real names of the probes (especially for the AHE data), while conforming to expected input of phyluce.

## Scripts used to convert


### find upper limit
[natsort](https://anaconda.org/anaconda/natsort) is not a native linux command
```
grep ">" Adephaga_2.9Kv1_UCE-Probes.fasta | cut -d"|" -f1 | cut -d_ -f1 | uniq | natsort | tail
```

Adephaga uce probes: highest ID **276883** 
Coleoptera no_intersect uce probes: 1172; 31213 total probes, **start at 500000**
vasil ahe probes: 49786 total probes, **start at 1000000**

### append real name to metadata
should be at end of header line and formated like this : `...,original-probe-ID:probeID` 
Note, I am unsure if the space between the start of the probe and pipe are necessary, BUT they will be kept to follow current formating standards
```
awk 'BEGIN{RS=">";-F" |"} (NR>1) {print ">" $1 " " $2",original-probe-ID:"$1 "\n" $3}' data/Coleoptera-UCE-1.1k_nointersect.fasta > Coleoptera-UCE-1.1K-v1_nointersect.tmp.fasta
awk 'BEGIN{RS=">";-F" |"} (NR>1) {print ">" $1 " " $2",original-probe-ID:"$1 "\n" $3}' data/vasilikopoulos_etal_ahe_nointersect.fasta > vasilikopoulos_etal_ahe_nointersect.tmp.fasta 
```

### replace probe with number
```
for i in 500000..500005..1; do awk 'BEGIN {-F" |"} (NR>1) {$1 = ">uce-'"$1"'_p$$??; print}' Coleoptera-UCE-1.1K-v1_nointersect.tmp.fasta; done | head

for i in {1..4..1}; do 
awk 'BEGIN {RS=">"; -F" |"} (NR>1) {$1 = ">uce-'"$i"'_p###"; print $1 $2 "\n" $3}' Coleoptera-UCE-1.1K-v1_nointersect.tmp.fasta; 
done | head
```


If I cannot get phyluce to use the data I am giving it we need to use a different pipeline.

https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html

blast, then use samtools
determine 

think about doing it by genome

ALSO LOOK AT A PAIRWISE JACARD DISTANCE BETWEEN ALL PROBES. For each genome e.g., 
*-GCA_911728475.2_icPteMadi1.2_genomic.fna.mem.sorted.mapped.bed
uce1 vs uce1, uce1 vs uce2, uce1 vs ahe, uce2 vs ahe, ahe vs ahe; see bedtools cookbook.


## probe maps for later concatenation

Need to generate a probe map for probes against a genome. Plan on using groupby to describe the intersection of probes on genes. I should be able to get a format that says what genomic feature has what probes:

Need to use the $PROBE-$GEN.mem.sorted.mapped.bed files against genome gff file

E.g., these files:
Adephaga_2.9Kv1_UCE-Probes.fasta-GCA_943142095.1_icOphArdo1.1_genomic.fna.mem.sorted.mapped.bed
Coleoptera-UCE-1.1K-v1.fasta-GCA_943142095.1_icOphArdo1.1_genomic.fna.mem.sorted.mapped.bed
vasilikopoulos_etal_ahe_probes.fasta-GCA_943142095.1_icOphArdo1.1_genomic.fna.mem.sorted.mapped.bed

against this genome gene feature format
GCA_943142095.1.gff

this should work similarly to the script found here:
https://github.com/crcardenas/Adephaga_UCE/blob/main/pterMadi2_script.sh

```
# sort GFF file
bedtools sort -i GCA_911728475.2.gff -g GCA_911728475.2_icPteMadi1.2_genomic.fna.genomefile > GCA_911728475.2.sorted.gff

# generate intergenic regions
awk '$3 == "gene" {print $0}' GCA_911728475.2.sorted.gff > GCA_911728475.2.sorted.genes.gff

# excise intergenic regions from gene features
bedtools complement -i GCA_911728475.2.sorted.genes.gff -g GCA_911728475.2_icPteMadi1.2_genomic.fna.genomefile > GCA_911728475.2.sorted.intergenic.bed

# extract exons
awk '$3 == "exon" {print $1 "\t" $4 "\t" $5}' GCA_911728475.2.sorted.gff > GCA_911728475.2.sorted.exons.bed # this will be so we can map intergenic regions
awk '$3 == "exon" {print $0}' GCA_911728475.2.sorted.gff > GCA_911728475.2.sorted.exons.gff # this ensures that we have a gff file to give intersect file exon names for later workflow

# join and sort
bedtools sort -i <(cat GCA_911728475.2.sorted.intergenic.bed GCA_911728475.2.sorted.exons.bed) -g GCA_911728475.2_icPteMadi1.2_genomic.fna.genomefile > GCA_911728475.2.intergenic-and-exons.sorted.bed

# excise intronic features 
bedtools complement -i GCA_911728475.2.intergenic-and-exons.sorted.bed -g GCA_911728475.2_icPteMadi1.2_genomic.fna.genomefile > GCA_911728475.2.sorted.introns.bed

```

###  Mapping probe to genome
	I am keeping ALL probes because I am not as concerned about specificity. I want to find what ever lands where regardless of if it was designed from this species. BWA is already specific enough (as compared to BLAST)
```
# index genome
bwa index GCA_911728475.2.fa  

# align probes to indexed genome
bwa mem GCA_911728475.2.fa Adephaga_2.9Kv1_UCE-Probes.fasta > probes.sam

# convert sam to bam
samtools view -h -b -S probes.sam > probes.mem.bam

# sort bam file
samtools sort probes.mem.bam -o probes.mem.sorted.bam

# extract sequences that map against the genome
samtools view -b -F 4 probes.mem.sorted.bam > probes.mem.sorted.mapped.bam

# create bed file
bedtools bamtobed -i probes.mem.sorted.mapped.bam > probes.mem.sorted.mapped.bed
```

# CHECK AGAIN HERE!

Next an intersect file is made
```
bedtools intersect -a GCA_911728475.2.gff \
    -b probes.mem.sorted.mapped.bed \
    -names cat4.4k \
#	-F 0.50 \
    -wa -wb > test.intersect
```

here we group the GFF feature info and sumarize by UCE loci
```
# -g group by column or ranges (1-4; 1,2,3,4)
# -c summarize by this column
bedtools groupby -i test.intersect -g 9 -c 13 -o collapse > test.intersect.grpby_gene
```
This identifies 1731 genes (not counting intergenic-genic)

You can now use a combination of grep/awk/cut and pipes to extract all the same information to used in [this R script](https://github.com/crcardenas/Adephaga_UCE/blob/main/just_pterostichus_probes.Rmd)
