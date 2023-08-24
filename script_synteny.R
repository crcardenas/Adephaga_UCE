
setwd("~/Desktop/CALOSOMA/SYNTENY/")

# Read in the Genome fasta index files

# Genome fasta index files
fai_sali<-read.table("GCA_944039245.1_icNebSali1.1_genomic.fna.fai")
fai_brev<-read.table("GCA_944738965.1_icNebBrev1.1_genomic.fna.fai")

# Add rownames for ordering the chromosomes:
rownames(fai_sali)<-fai_sali$V1
rownames(fai_brev)<-fai_brev$V1

# Remove tiny scaffolds
#fai_sali<-fai_sali[!grepl(fai_sali$V1,pattern="unloc|scaffold|:"),]
#fai_brev<-fai_brev[!grepl(fai_brev$V1,pattern="unloc|scaffold|:"),]

# Order according to chromosome names
require(gtools)
fai_sali<-fai_sali[mixedsort(rownames(fai_sali),decreasing = F),]
fai_brev<-fai_brev[mixedsort(rownames(fai_brev),decreasing = F),]

# # Order according to M. menophilus chromosome numbers
#fai_meno<-fai_meno[c("SUPER_9","SUPER_16","SUPER_12","SUPER_10","SUPER_5","SUPER_4","SUPER_13","SUPER_7","SUPER_14","SUPER_3","SUPER_8","SUPER_17","SUPER_15","SUPER_19","SUPER_1","SUPER_18","SUPER_11","SUPER_2","SUPER_6","SUPER_20","SUPER_Z"),]

# Add additive positions (to plot chromosomes next to each other)
fai_sali$add<-cumsum(c(0,fai_sali$V2[-length(fai_sali$V2)]+11000000))
fai_brev$add<-cumsum(c(0,fai_brev$V2[-length(fai_brev$V2)]+7700000))

# Add a colour column where each chromosome has a different colour using a colourblind friendly palette
require(colorBlindness)
fai_sali$col<-paletteMartin[1:nrow(fai_sali)]

# Read in BUSCO gene coordinates
busco_brev<-read.table("C_NebBrev.bed")
busco_sali<-read.table("C_NebSali.bed")

# Rename the columns
names(busco_brev)<-c("chr_brev","start_brev","end_brev","gene")
names(busco_sali)<-c("chr_sali","start_sali","end_sali","gene")

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
