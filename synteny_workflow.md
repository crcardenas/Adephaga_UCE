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
`CHECK FAI FILES, one of the species loaded into R while testing has a strange count?`

### Data is ready for R script `script_synteny.R`

CHECK FAI FILES

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
## Modified from ? Jérémy Gauthier's synteny_script.R

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
  print(paste0("Genome input:  ",))
  print("")
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
                           paste0("loci_",species1))

sp2_bedfile <- read.table(c(paste0(species2,".ortho-loci.bed")))
colnames(sp2_bedfile) <- c(paste0("chr_",species2),
                           paste0("start_",species2),
                           paste0("end_",species2),
                           paste0("loci_",species2))

# genome fai index files
sp1_fai <- read.table(c(paste0("../",species1,".fna.fai")))
rownames(sp1_fai) <- sp1_fai$V1

sp2_fai <- read.table(c(paste0("../",species2,".fna.fai")))
rownames(sp2_fai) <- sp2_fai$V1

#### Manipulate data ####
# Order according to chromosome names
sp1_fai <- sp1_fai[mixedsort(rownames(sp1_fai),decreasing = F),]
sp2_fai <- sp2_fai[mixedsort(rownames(sp2_fai),decreasing = F),]

#### scripting notes ####
# use some variable in file path
# dir <- file.path("some", "path")
# bla <- file.path("some", "directory")
# files <- c("file1.R", "file2.exe")
# 
# file.path(dir, bla, files)
#
# results in:
# [1] "some/path/some/directory/file1.R"   "some/path/some/directory/file2.exe"
# 
# might be able to call paste0(c()) within a function. instead of:
# workingdir <- paste0("./",species1,"_",species2,"/")
# *maybe*
# myfunciton(directory=paste0("./",species1,"_",species2,"/"), other_arg=0, ...)


#### Jeremy's code ####
# !!! if commented out that means I've added that to my modified script


   # Read in the Genome fasta index files
  
   # Genome fasta index files
   # fai_sali<-read.table("GCA_944039245.1_icNebSali1.1_genomic.fna.fai")
   # fai_brev<-read.table("GCA_944738965.1_icNebBrev1.1_genomic.fna.fai")
  
   # Add rownames for ordering the chromosomes:
   # rownames(fai_sali)<-fai_sali$V1
   # rownames(fai_brev)<-fai_brev$V1
  
   # Remove tiny scaffolds
   #fai_sali<-fai_sali[!grepl(fai_sali$V1,pattern="unloc|scaffold|:"),]
   #fai_brev<-fai_brev[!grepl(fai_brev$V1,pattern="unloc|scaffold|:"),]
  
   # # Order according to chromosome names
   # require(gtools)
   # fai_sali<-fai_sali[mixedsort(rownames(fai_sali),decreasing = F),]
   # fai_brev<-fai_brev[mixedsort(rownames(fai_brev),decreasing = F),]
  
   # # Order according to M. menophilus chromosome numbers
   #fai_meno<-fai_meno[c("SUPER_9","SUPER_16","SUPER_12","SUPER_10","SUPER_5","SUPER_4","SUPER_13","SUPER_7","SUPER_14","SUPER_3","SUPER_8","SUPER_17","SUPER_15","SUPER_19","SUPER_1","SUPER_18","SUPER_11","SUPER_2","SUPER_6","SUPER_20","SUPER_Z"),]
  
  # OK porquoi? why add length?
   # Add additive positions (to plot chromosomes next to each other)
   fai_sali$add<-cumsum(c(0,fai_sali$V2[-length(fai_sali$V2)]+11000000))
   fai_brev$add<-cumsum(c(0,fai_brev$V2[-length(fai_brev$V2)]+7700000))
  
  # OK porquoi? for... only one??
   # Add a colour column where each chromosome has a different colour using a colourblind friendly palette
   require(colorBlindness)
   fai_sali$col<-paletteMartin[1:nrow(fai_sali)]
  
   # # Read in BUSCO gene coordinates
   # busco_brev<-read.table("C_NebBrev.bed")
   # busco_sali<-read.table("C_NebSali.bed")
   # 
   # # Rename the columns
   # names(busco_brev)<-c("chr_brev","start_brev","end_brev","gene")
   # names(busco_sali)<-c("chr_sali","start_sali","end_sali","gene")
  
   # remove tiny scaffolds from the BUSCO gene set
   #busco_mars<-busco_mars[!grepl(busco_mars$chr_mars,pattern="unloc|scaffold|:"),]
   #busco_meno<-busco_meno[!grepl(busco_meno$chr_meno,pattern="unloc|scaffold|:"),]
  
   # swap orientation of some M. menophilus chromosomes to better match the M. marsaeus chr
   chrLengths<-as.data.frame(aggregate(busco_brev$end_brev ~ busco_brev$chr_brev,FUN = max))
   names(chrLengths)<-c("chr_brev","length")
   for(chr in chrLengths$chr_brev){
     busco_brev[busco_brev$chr_brev==chr,2:3]<-chrLengths[chrLengths$chr_brev==chr,2]-busco_brev[busco_brev$chr_brev==chr,2:3]
   }
   
   # Add a column with additive positions
   busco_sali$addPos_sali<-rowMeans(cbind(busco_sali$start_sali,busco_sali$end_sali))+fai_sali[busco_sali$chr_sali,"add"]
   busco_brev$addPos_brev<-rowMeans(cbind(busco_brev$start_brev,busco_brev$end_brev))+fai_brev[busco_brev$chr_brev,"add"]
  
  
   # Merge the two gene coordinate tables
   busco<-merge(busco_sali,busco_brev,by="gene")
   busco2<-busco[order(busco$chr_sali,busco$start_sali),]
  
   # Write out and get the correct order in Excel
   # write.table(busco,file="D:/Dropbox/Ithomiines/ReferenceGenomes/curatedMelinaeaGenomes/busco.comparison.txt")
   # busco<-read.table("D:/Dropbox/Ithomiines/ReferenceGenomes/curatedMelinaeaGenomes/busco.comparison.txt",header=T)
  
   # Plot the gene coordinates against each other
   par(mfrow=c(1,1),mar=c(4,4,1,1),mgp=c(1.7,0.5,0),xaxs="i",yaxs="i")
   plot(busco2$addPos_sali,busco2$addPos_brev,xaxt="n",yaxt="n",
        col=as.integer(as.factor(busco2$chr_sali)),
        pch=19,cex=0.5,
        xlab=expression(italic("Nebria brev")~"chromosomes"),
        ylab=expression(italic("Nebria sali")~"chromosomes"))
  
   plot(busco$addPos_sali,busco$addPos_brev,xaxt="n",yaxt="n",
        col=as.integer(as.factor(busco$chr_sali)),
        pch=19,cex=0.5,xlab=expression(italic("Nebria brev")~"chromosomes"),
        ylab=expression(italic("Nebria sali")~"chromosomes"))
  
  
  
   # Add lines showing the ends of the chromosomes
   abline(v=fai_sali$add[-1][!grepl(fai_sali$V1,pattern="unloc")],col="grey")
   abline(h=fai_brev$add[-1],col="grey")
  
   # Label the axes
   axis(1,fai_sali$add[!grepl(fai_sali$V1,pattern="unloc")]+
          fai_sali$V2[!grepl(fai_sali$V1,pattern="unloc")]/2,
        fai_sali$V6[!grepl(fai_sali$V1,pattern="unloc")],line = -0.2,tick = F)
   axis(2,fai_brev$add+fai_brev$V2/2,fai_brev$V6,las=2,tick = F)


```
