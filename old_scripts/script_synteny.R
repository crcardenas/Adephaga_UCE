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
exclude <- unique(c(sp1_exclude,sp2_exclude))
# remove those loci on those scaffolds
sp1_bedfile <- sp1_bedfile[ ! (sp1_bedfile[,4] %in% exclude) , ]
sp2_bedfile <- sp2_bedfile[ ! (sp2_bedfile[,4] %in% exclude) , ]

# Order according to chromosome names
sp1_fai <- sp1_fai[mixedsort(rownames(sp1_fai),decreasing = F),]
sp2_fai <- sp2_fai[mixedsort(rownames(sp2_fai),decreasing = F),]

# add space between the chromosomes for visualisation
sp1_fai$add<-cumsum(c(0,sp1_fai$V2[-length(sp1_fai$V2)]+11000000))
sp2_fai$add<-cumsum(c(0,sp2_fai$V2[-length(sp2_fai$V2)]+7000000))

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
        xlab=c(paste0(species2," chromosomes")),
        ylab=c(paste0(species1," chromosomes")))
  
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
