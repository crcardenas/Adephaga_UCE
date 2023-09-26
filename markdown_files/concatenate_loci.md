# Concatenate Loci for phylogenetic analysis

Two choices, concatenate loci that are found on the same gene OR concatenate loci based on their distance from one another. The assumption for by gene is that loci on the gene should have a similar evolutionary rate. However, genes can have varying rates of evolution within them. Creating a partition by length, assuming linkage between loci, it is possible to join loci that need joined. This, is mainly for demonstration.


1) map and make bedfile 
2) find clusters using `bedtools cluster -d` and check a range from 0:1000 and examine the change in the number of "clusters"
3) write script to concatenate loci/uce's by the cluster statistic chosen



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
	GENE=$(echo $TMP | cut -d ";" -f 1 | cut -d ":" -f 2)
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
this gives a filethat looks like 

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

First make a directory to work from `7modify-partition`, and copy your phyluce partition file into it. Because phyluce uses the file name, and not just the loci/uce name, it needs the `*.nexus` removed from the file. Additionally, we only want to rename UCEs if there is more than one found on a gene. Depending on the completeness of the matrix that gets put in, then many taxa might be missing one or more of this loci/probe 

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
The new partition can be run in phyluce