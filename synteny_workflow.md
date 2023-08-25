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

`CLEAN UP THIS CODE, REDUCE REDUNDANCY`
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

`CLEAN UP THIS CODE, REDUCE REDUNDANCY`
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

`CLEAN UP THIS CODE, REDUCE REDUNDANCY`
`ELABORATE ON WHAT THIS STEP DOES`

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

`Bash script should work`

### `synteny.sh`
```
#!/bin/bash

species=("nebBrev1" "nebSali1" "ophArdo1" "pteMadi2")

for ((i=0; i<${#species[@]}; i++)); do
    for ((j=i+1; j<${#species[@]}; j++)); do

		species_1="${species[i]}"
    	species_2="${species[j]}

		# preprocessing steps

		# run R.r script with ${_output_file} parameter
		Rscript --vanilla script_synteny.R "${species_1}" "${species_2}"

		#printf "synteny plot created for loci of $species_1 vs $species_2""

		exit 0
	done
done
```

### `script_synteny.R`

```
#### script_synteny.R ####
## Cody Raul Cardenas 
## 2023.08 
## Modified from Joana Meier's synteny_script.R

# to be run with a bash script like this
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~                           synteny.sh                          ~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #!/bin/bash
# 
# species=("nebBrev1" "nebSali1" "ophArdo1" "pteMadi2")
# 
# for ((i=0; i<${#species[@]}; i++)); do
#   for ((j=i+1; j<${#species[@]}; j++)); do
#     
#     species_1="${species[i]}"
#     species_2="${species[j]}"
#     
#     # preprocessing steps
#     
#     Rscript --vanilla script_synteny.R "${species_1}" "${species_2}"
#     
#     done
# done
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### house keeping script ####
# get the input passed from the shell script
args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Two argument must be supplied as input.\n", call. = FALSE)
} else {
  print(paste0("Genome input:  ", args[1]))
  print(paste0("Genome input:  ", args[2]))
  print("Saving file")
  print("Process output in excel or equivelant")
}

# load library
library(gtools)
library(colorBlindness)

#### user input scripts ####
# set variables load data
# species
species1 <- args[1]
species2 <- args[2]

# working directory
setwd(c(paste0(species1,"_",species2)))

# bedfiles
sp1_bedfile <- read.table(c(paste0(species1,".ortho-loci.bed")))
colnames(sp1_bedfile) <- c(paste0("chr_",species1),
                           paste0("start_",species1),
                           paste0("end_",species1),
                           paste0("loci"))

sp2_bedfile <- read.table(c(paste0(species2,".ortho-loci.bed")))
colnames(sp2_bedfile) <- c(paste0("chr_",species2),
                           paste0("start_",species2),
                           paste0("end_",species2),
                           paste0("loci"))

# genome fai index files
sp1_fai <- read.table(c(paste0("../",species1,".fna.fai")))
rownames(sp1_fai) <- sp1_fai$V1

sp2_fai <- read.table(c(paste0("../",species2,".fna.fai")))
rownames(sp2_fai) <- sp2_fai$V1

#### Manipulate data ####
# retain only chromosomes
sp1_fai <- subset(sp1_fai, nchar(as.character(V1)) <= 3)
sp2_fai <- subset(sp2_fai, nchar(as.character(V1)) <= 3)

# remove any UCE loci from these chromosomes
# make a list to check if loci on scaffolds
sp1_exclude <- subset(sp1_bedfile, nchar(sp1_bedfile[,1]) > 3)[,4]
sp2_exclude <- subset(sp2_bedfile, nchar(sp2_bedfile[,1]) > 3)[,4]
# create a list
exclude <- paste0(sp1_exclude,sp2_exclude)
# remove those loci on those scaffolds
sp1_bedfile <- sp1_bedfile[ ! (sp1_bedfile[,4] %in% exclude) , ]
sp2_bedfile <- sp2_bedfile[ ! (sp2_bedfile[,4] %in% exclude) , ]

# Order according to chromosome names
sp1_fai <- sp1_fai[mixedsort(rownames(sp1_fai),decreasing = F),]
sp2_fai <- sp2_fai[mixedsort(rownames(sp2_fai),decreasing = F),]

# add space between the chromosomes for visualisation
sp1_fai$add<-cumsum(c(0,sp1_fai$V2[-length(sp1_fai$V2)]+11000000))
sp2_fai$add<-cumsum(c(0,sp2_fai$V2[-length(sp2_fai$V2)]+7700000))

# add color for plotting
sp1_fai$col<-paletteMartin[1:nrow(sp1_fai)]

# WORK ON THESE CODES
# FOR TESTING swap orientation of chromosomes in species one
chrLengths_sp1 <- as.data.frame(aggregate(sp1_bedfile[,3] ~ sp1_bedfile[,1],
                                      FUN = max))
colnames(chrLengths_sp1) <- c(paste0("chr_",species1),
                           paste0("len_",species1))
for ( chr in chrLengths_sp1[,1] ) {
  sp1_bedfile[sp1_bedfile[,1]==chr,2:3] <- chrLengths_sp1[chrLengths_sp1[,1]==chr,2] - sp1_bedfile[sp1_bedfile[,1]==chr,2:3]
}

# Add a column with additive positions
sp1_bedfile$V5 <- rowMeans(cbind(sp1_bedfile[,2], sp1_bedfile[,3])) + sp1_fai[sp1_bedfile[,1], "add"]
colnames(sp1_bedfile) <- c(paste0("chr_",species1),
                          paste0("start_",species1),
                          paste0("end_",species1),
                          paste0("loci"),
                          paste0("addPos_",species1))
  
sp2_bedfile$V5 <- rowMeans(cbind(sp2_bedfile[,2], sp2_bedfile[,3])) + sp2_fai[sp2_bedfile[,1], "add"]
colnames(sp2_bedfile) <- c(paste0("chr_",species2),
                           paste0("start_",species2),
                           paste0("end_",species2),
                           paste0("loci"),
                           paste0("addPos_",species2))

# Merge the two gene coordinate tables
loci <- merge(sp1_bedfile, sp2_bedfile, by="loci")
loci2 <- loci[order(loci[,6], loci[,7]),]

# MAY NEED TO WRITE TWO SCRIPTS
# ONE TO GET THIS FILE (loci1, loci2)
# THEN THAT FILE NEEDS TO BE MANUALLY CURRATED IN EXCEL or EQUIVELANT

#    Write out and get the correct order in Excel
#    write.table(busco,file="D:/Dropbox/Ithomiines/ReferenceGenomes/curatedMelinaeaGenomes/busco.comparison.txt")
#    busco<-read.table("D:/Dropbox/Ithomiines/ReferenceGenomes/curatedMelinaeaGenomes/busco.comparison.txt",header=T)


#### plot ####
#pdf(file=c(paste0(species1,"_",species2,"_syntenyplot.pdf")))
   # Plot the gene coordinates against each other
   par(mfrow=c(1,1),mar=c(4,4,1,1),mgp=c(1.7,0.5,0),xaxs="i",yaxs="i")
   plot(loci2[,5],
        loci2[,9],
        xaxt="n",
        yaxt="n",
        col=as.integer(as.factor(loci2[,6])),
        pch=19,cex=0.5,
        xlab=c(paste0(species1," chromosomes")),
        ylab=c(paste0(species2," chromosomes")))
  
   # plot(loci[,9],
   #      loci[,5],
   #      xaxt="n",yaxt="n",
   #      col=as.integer(as.factor(loci2[,6])),
   #      pch=19, cex=0.5, 
   #      xlab=c(paste0(species1," chromosomes")),
   #      ylab=c(paste0(species2," chromosomes")))
   # 
  
   # Add lines showing the ends of the chromosomes
   abline(v=sp1_fai$add[-1][!grepl(sp1_fai$V1,pattern="unloc")],col="grey")
   abline(h=sp2_fai$add[-1],col="grey")

  # Label the axes
   # axis(2, sp1_fai[,1], line = -0.2, tick = F)
   # axis(1, sp2_fai[,1], las = 2, tick = F)
   
# not sure what this did...
   # axis(1,sp2_fai$add[!grepl(sp2_fai$V1, pattern="unloc")]+
   #        sp2_fai$V2[!grepl(sp2_fai$V1, pattern="unloc")]/2,
   #      sp2_fai$V6[!grepl(sp2_fai$V1, pattern="unloc")], line = -0.2, tick = F)
   # axis(2, sp1_fai$add+sp1_fai$V2/2, sp1_fai$V6,las=2,tick = F)
#dev.off()

```
