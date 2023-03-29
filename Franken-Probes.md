---
title: Franken Probes
author: CR Cardenas
date: 2023 03 29
---

# Generate non-intersecting probes

Need to map multiple genomes to realize the overalp of all probesets. Remember AHE probes are like this TC007324-PA_1 is the gene TC007324-PA_1_4 is the exon, all other codes are probe identifiers

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

code requires the following software: samtools, bedtools, and bwa. 
useful tools used but not necessary: seqkit

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
PROBE=$2

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


