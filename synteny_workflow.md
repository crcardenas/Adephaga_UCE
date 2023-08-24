---
title: Synteny Plotting
author: CR Cardenas
date: 2023.08
---
# Synteny Plotting of UCEs
goal is to identify how UCEs have shifted around genomes.
1. find shared UCEs parwise by species
	1. e.g., sp1 x sp2, sp1 x sp3, sp1 x sp4, sp2 x sp3, sp2 x sp4, sp3 x sp4
2. create a samtools index of genomes
3. use R script provided by Jeremy to create a synteny plot for each species
## input data
* a genome fasta file/s
* a intersect file/s from adephaga workflow
* probe file/s
# Data processing
First create a bed file of just the chromosome position, the start, the end, and probe name from the intersect files in dir `3intersect`
```
# create simplified bedtools file
awk -F"\t" '{print $11"\t"$12"\t"$13"\t"$14}' pteMadi2-PROBES.intersect > pteMadi2-PROBES.intersect.bed

#sort
bedtools sort -i pteMadi2-PROBES.intersect.bed > pteMadi2-PROBES.intersect.sorted.bed

awk -F"\t" '{print $11"\t"$12"\t"$13"\t"$14"\t"$10}' ophArdo1-PROBES.intersect | bedtools sort > ophArdo1-PROBES.intersect.sorted.bed
```
### rename Vasil and Coleoptera UCEs (see workflow for Adephaga UCE characterization)
so we that all of the loci being used for the phylogenetic anlaysis can be examined for synteny, we need to rename the loci so they are easy to parse and dont get confused between probesets. Ensure that the Unaltered probe files are in the same directory. I am using ALL probes rather than the subsets!!!

create a list based on fasta
```
grep ">" Coleoptera.fasta | cut -d"|" -f1 | tr -d ">" > Coleoptera-rename.list

grep ">" Vasil.fasta | cut -d"|" -f1 | tr -d ">" > Vasil-rename.list

touch Coleoptera-rename.tmp.list
for i in $(cat Coleoptera-rename.list); do
LOCI=$(echo $i | cut -d "p" -f 1);
echo -e $i "\t" $LOCI >> Coleoptera-rename.tmp.list;
done

touch Vasil-rename.tmp.list
for i in $(cat Vasil-rename.list); do
LOCI=$(echo $i | cut -d _ -f 1,2,3);
echo -e $i "\t" $LOCI >> Vasil-rename.tmp.list; 
done

awk 'BEGIN {count=1000000;probe[0]=""};
# if the probe is not loaded into the array; that probe gets +1
{if (!($2 in probe)) {probe[$2]=count++} ;
# print a new column with your new loci name
print $0, "\t", "uce-"count}' Coleoptera-rename.tmp.list > Coleoptera-rename.tmp.list2

# do the same, but for our newly created column
# compare the value from the first field stored in save
awk '$3 != save { counter = 1; save = $3 }
# if the values differ resets to one
{ print $0, "\t", "_p"counter++ }' Coleoptera-rename.tmp.list2 > Coleoptera-rename.tmp.list3

# and kiss
awk '{print $1, "\t", $3 $4}' Coleoptera-rename.tmp.list3 > Coleoptera.renamed-probe.list

#again for Vasil but at 2 million so there is no overlap
awk 'BEGIN {count=2000000;probe[0]=""};
{if(!($2 in probe)) {probe[$2]=count++};
print $0, "\t", "uce-"count}' Vasil-rename.tmp.list > Vasil-rename.tmp.list2
awk '$3 != save { counter = 1; save = $3 } { print $0, "\t", "_p" counter++ }' Vasil-rename.tmp.list2 > Vasil-rename.tmp.list3
awk '{print $1, "\t", $3 $4}' Vasil-rename.tmp.list3 > Vasil.renamed-probe.list
```
substitute your bedfiles names for both the coleoptera and vasil files
```
# substitute for coleoptera names
awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Coleoptera") { $4 = map[$4] } {print $0}' Coleoptera.renamed-probe.list pteMadi2-PROBES.intersect.sorted.bed > pteMadi2.tmp
# substitute for vasil names
awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Vasil") { $4 = map[$4] } {print $0}' Vasil.renamed-probe.list pteMadi2.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4}' > pteMadi2.renamed.bed

awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Coleoptera") { $4 = map[$4] } {print $0}' Coleoptera.renamed-probe.list nebBrev1-PROBES.intersect.sorted.bed > nebBrev1.tmp
awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Vasil") { $4 = map[$4] } {print $0}' Vasil.renamed-probe.list nebBrev1.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4}' > nebBrev1.renamed.bed

awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Coleoptera") { $4 = map[$4] } {print $0}' Coleoptera.renamed-probe.list nebSali1-PROBES.intersect.sorted.bed > nebSali1.tmp
awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Vasil") { $4 = map[$4] } {print $0}' Vasil.renamed-probe.list nebSali1.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4}' > nebSali1.renamed.bed

awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Coleoptera") { $4 = map[$4] } {print $0}' Coleoptera.renamed-probe.list ophArdo1-PROBES.intersect.sorted.bed > ophArdo1.tmp
awk 'FNR==NR { map[$1] = $2; next } $4 in map && ($5="Vasil") { $4 = map[$4] } {print $0}' Vasil.renamed-probe.list ophArdo1.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ophArdo1.renamed.bed
```
## Find start and end of individual Loci
### `loci_positions.awk`
To run
```
awk -f loci_positions.awk nebBrev1.renamed.bed | bedtools sort -i > nebBrev1-LOCI.bed
awk -f loci_positions.awk nebSali1.renamed.bed | bedtools sort -i > nebSali1-LOCI.bed
awk -f loci_positions.awk pteMadi2.renamed.bed | bedtools sort -i > pteMadi2-LOCI.bed
awk -f loci_positions.awk ophArdo1.renamed.bed | bedtools sort -i > ophArdo1-LOCI.bed
```
the code
### `CHECK THE CODE HERE, IT INCLUDES SOME PROBE DATA IN THE FINAL CODE, CONSIDER CHANGING THE SPLIT TO 'p' SO THE STRING IS CONSISTENT AND THERE IS NO CONFUSING uce-123 WITH uce-1234`
`I've manually fixed the output bed files..., but it still needs cleaned up`
```
#!/usr/bin/awk -f
# The BEGIN block sets the field separator (FS) and output field separator (OFS) for processing tab-separated values.
BEGIN {FS="\t"; OFS="\t"}

# For each line in the input, split the fourth column (loci information) using underscore as the separator.
{
    split($4, loci_parts, "_");
    loci_id = loci_parts[1];  # Extract the loci ID without the probe ID.
    
    # If the loci ID is not already in the loci_data associative array, store the entire line.
    if (!(loci_id in loci_data)) {
        loci_data[loci_id] = $0;
    } else {
        # If the loci ID already exists in loci_data, split the stored line and compare/merge the start and end positions.
        split(loci_data[loci_id], prev, "\t");
        prev[2] = (prev[2] < $2) ? prev[2] : $2;  # Compare and update the start position.
        prev[3] = (prev[3] > $3) ? prev[3] : $3;  # Compare and update the end position.
        loci_data[loci_id] = prev[1] OFS prev[2] OFS prev[3] OFS loci_id;  # Reconstruct the line with merged positions.
    }
}

# After processing all lines, loop through the loci_data associative array and print the merged lines.
END {
    for (loci_id in loci_data) {
        print loci_data[loci_id];
    }
}
```
# Identify common loci between bed files
First identify the common loci between bedfiles
second create a new file of "orthologous" loci (at least shared between these species)
### `common_loci.sh`
```
#!/bin/bash

species=("nebBrev1" "nebSali1" "ophArdo1" "pteMadi2")

for ((i=0; i<${#species[@]}; i++)); do
    for ((j=i+1; j<${#species[@]}; j++)); do
        
		species_1="${species[i]}"
        species_2="${species[j]}"

		newdirectory=${species_1}_${species_2}

		output_file="common_${species_1}_${species_2}.list"

		if [ ! -d ${newdirectory} ]; then

			mkdir ${newdirectory}

			awk -v species_1="$species_1" -v species_2="$species_2" 'BEGIN {FS="\t"; OFS="\t"} {
				loci_name = $4;
				if (FILENAME == ARGV[1]) {
					loci_sp1[loci_name] = 1;
				} else if (FILENAME == ARGV[2]) {
					if (loci_name in loci_sp1) {
						print $4;
					}
				}
			}' "${species_1}-LOCI.bed" "${species_2}-LOCI.bed" > "$newdirectory"/"$output_file"
			
			cp "${species_1}-LOCI.bed" "$newdirectory"/
			cp "${species_2}-LOCI.bed" "$newdirectory"/
		

        else 

			awk -v species_1="$species_1" -v species_2="$species_2" 'BEGIN {FS="\t"; OFS="\t"} {
				loci_name = $4;
				if (FILENAME == ARGV[1]) {
					loci_sp1[loci_name] = 1;
				} else if (FILENAME == ARGV[2]) {
					if (loci_name in loci_sp1) {
						print $4;
					}
				}
			}' "${species_1}-LOCI.bed" "${species_2}-LOCI.bed" > "$newdirectory"/"$output_file"
			
			cp "${species_1}-LOCI.bed" "$newdirectory"/
			cp "${species_2}-LOCI.bed" "$newdirectory"/

		fi

    done

done

for ((i=0; i<${#species[@]}; i++)); do
	for ((j=i+1; j<${#species[@]}; j++)); do
	
	species_1="${species[i]}"
    species_2="${species[j]}"
	
	directory="${species_1}_${species_2}"

	sp1_bed="${species_1}-LOCI.bed"
	sp2_bed="${species_2}-LOCI.bed"
	
	list="common_${species_1}_${species_2}.list"
	
	sp1_output="${species_1}.ortho-loci.bed"
	sp2_output="${species_2}.ortho-loci.bed"

	awk 'FNR==NR{bed[$1]; next} $4 in bed' "${directory}/${list}" "${directory}/${sp1_bed}" > "${directory}/${sp1_output}"
	awk 'FNR==NR{bed[$1]; next} $4 in bed' "${directory}/${list}" "${directory}/${sp2_bed}" > "${directory}/${sp2_output}"
	
	echo "Pairwise comparison between $species_1 and $species_2 completed. Output saved to $directory/$sp1_output & $directory/$sp2_output"

	done

done
```
# Create samtools index files
```
for i in *.fna; do samtools faidx $i; done
``` 

### Data is ready for R script `script_synteny.R`

