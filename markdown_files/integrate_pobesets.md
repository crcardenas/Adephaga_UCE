# Integrate probesets

If you find it necessary to join multiple probesets follow the rest of this pipeline. Some steps of this pipeline will provide a summary of the probes that may be useful.

## Data clean up

Given that Adephaga has ~300 UCE loci in common with Coleoptera by design, these will be removed from the Adephaga probekit. This pipeline is prioritizing Adephaga, but either probekit could be used to remove. Given that the original Coleoptera probekit is being integrated, those loci/probes will be removed from the Adephaga set.

This reduces the 2921 targeted loci from the previous steps to 2619:
```
# get list of probes to remove
grep "coleoptera" 0data/probes/Adephaga_2.fasta | cut -d "|" -f 1 | cut -d "p" -f 1 | sort -u > 0data/probes/tmp.list

# use list to remove those probes
awk '{ if ((NR>1) && ($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' 0data/probes/Adephaga_2.fasta | grep -vf 0data/probes/tmp.list - | tr "\t" "\n" > 0data/probes/Adephaga_3.fasta
```

## Subset overlap between probefiles 

### `blast.sh`

order is important
if you have an probefile that isn't in the expected phyluce format (e.g., `uce100_p10 | ...`) then you need to modify the cut in the last awk statements or do this step manually. This file takes input [....]

needs to be run twice;
```
bash blast.sh Adephaga_3.fasta Coleoptera_2.fasta
bash blast.sh Adephaga_3.fasta Vasil.fasta
```


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
	cut -d "p" -f 1 | \
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
	cut -d "p" -f 1 | \
	sort -u \
	> '4probe-reduction/'$PROBEFILE2NAME'.to-remove'

fi
```

Next add to the `*.to-remove` list based on intersect files generated in directory `3intersect.` 

```
awk '$10 == "Coleoptera" {print $14}' 3intersect/pteMadi2-PROBES.intersect | cut -d "p" -f 1 | sort -u >> 4probe-reduction/Coleoptera_2.to-remove
awk '{print $0}' 4probe-reduction/Coleoptera_2.to-remove | sort -u > 4probe-reduction/Coleoptera_2.to-remove.list 
```

In this example, the Vasil probes are formated differently and for downstream purposes I am considering exons over the genes (note the change in cut)

```
awk '$7 == "Vasil" {print $11}' 3intersect/pteMadi2-PROBES.intersect | cut -d _ -f 1,2,3 | sort -u >> 4probe-reduction/Vasil.to-remove2
awk '{print $0}' 4probe-reduction/Vasil.to-remove2 | sort -u > 4probe-reduction/Vasil.to-remove.list
```

next we remove the duplicates using the generated lists

```
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


```
touch 5join-probes/Coleoptera-rename.tmp.list
for i in $(cat 5join-probes/Coleoptera-rename.list); do
LOCI=$(echo $i | cut -d "p" -f 1);
echo -e $i "\t" $LOCI >> 5join-probes/Coleoptera-rename.tmp.list;
done

touch 5join-probes/Vasil-rename.tmp.list
for i in $(cat 5join-probes/Vasil-rename.list); do
LOCI=$(echo $i | cut -d _ -f 1,2,3);
echo -e $i "\t" $LOCI >> 5join-probes/Vasil-rename.tmp.list; 
done

```

Next create a new ID for these probes based on their "loci" to make it compatable with `phyluce`. This uses the new tmp.list file created previously to alter the fasta file.

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
awk '$3 != save { counter = 1; save = $3 } { print $0, "\t", "_p" counter++ }' 5join-probes/Vasil-rename.tmp.list2 > 5join-probes/Vasil-rename.tmp.list3
awk '{print ">" $1, "\t", ">" $3 $4}' 5join-probes/Vasil-rename.tmp.list3 > 5join-probes/Vasil.renamed-probe.list

# clean up tmp files
#rm 5join-probes/*.tmp.list*
```

## Rename and concatenate probes

First we rename the probe files using the `*.renamed-probe.list` and the `*.subset.fasta`

```
awk 'BEGIN {OFS="\t"} NR==FNR {a[$1]=$2; next} /^>/ {if (a[$1]) {$1=a[$1]}}1' 5join-probes/Coleoptera.renamed-probe.list 5join-probes/Coleoptera.subset.tmp  | tr "\t" " " > 5join-probes/Coleoptera.renamed.tmp.probes

awk 'BEGIN {OFS="\t"} NR==FNR {a[$1]=$2; next} /^>/ {if (a[$1]) {$1=a[$1]}}1' 5join-probes/Vasil.renamed-probe.list 5join-probes/Vasil.subset.tmp  | tr "\t" " " > 5join-probes/Vasil.renamed.tmp.probes

cat 0data 0data/probes/Adephaga_3.fasta 5join-probes/Coleoptera.renamed.tmp.probes 5join-probes/Vasil.renamed.tmp.probes > 5join-probes/joined.tmp.fasta
```

following this pipeline the resulting probeset contains 4255 loci comprised of 90130 probes. So the new probefile will follow the naming convention of other UCE probesets joined_probes4.2k.fasta

```
mv 5join-probes/joined.renamed.tmp.probes > 5join-probes/joined_probes4.2k.fasta
```

Here are the resulting changes for the probes

`$ seqkit stats *.fasta`

|file|format|type|num_seqs|sum_len|min_len|avg_len|max_len|
|---|---|---|---|---|---|---|---|
|0data/probes/Adephaga.fasta|FASTA|DNA|38,948|4,673,760|120|120|120|
|0data/probes/Adephaga_2.fasta|FASTA|DNA|38,682|4,641,840|120|120|120|
|0data/probes/Adephaga_3.fasta|FASTA|DNA|34,995|4,199,400|120|120|120|
|0data/probes/Coleoptera.fasta|FASTA|DNA|13,674|1,640,880|120|120|120|
|0data/probes/Coleoptera_2.fasta|FASTA|DNA|13,609|1,633,080|120|120|120|
|0data/probes/Vasil.fasta|FASTA|DNA|49,786|5,974,320|120|120|120|
|5join-probes/Coleoptera.subset.tmp|FASTA|DNA|9,033|1,083,960|120|120|120|
|5join-probes/Vasil.subset.tmp|FASTA|DNA|46,102|5,532,240|120|120|120|
|5join-probes/joined_probes4.2k.fasta|FASTA|DNA|90,130|10,815,600|120|120|120|


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
