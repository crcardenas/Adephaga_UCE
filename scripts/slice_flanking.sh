#!/bin/bash
source /local/anaconda3/bin/activate
# this conda environment already has all the software & packages necessary
conda activate characterization

# test BWA index
mkdir tmp

cd tmp
for i in ../loci/*.fasta; do
echo indexing ${i}

# add indexing logic to skip indexing if completed (if ... else ... fi statement)
# index uce loci
	bwa index ${i}; 

echo -e "\n"	
done

for i in ../loci/*.fasta; do

FASTA=$(echo $i | cut -d "/" -f 3)
echo mapping ${FASTA}

# map ID probe region in loci, data from phyluce slice function using: -flank 0 (noflank) from -flank 400 (loci)
	bwa mem -t 6 ${i} ../noflank/${FASTA} \
		> ./${FASTA}.mapped.sam;
# generate genome file
	samtools faidx ${i};
# convert to bam file
	samtools view -h -b -S ./${FASTA}.mapped.sam \
		> ./${FASTA}.mapped.bam;
# convert to bed file
	bedtools bamtobed -i ./${FASTA}.mapped.bam \
		> ../bed_files/${FASTA}.mapped.bed;

echo -e "\n"
done

for i in ../loci/*.fasta; do
FASTA=$(echo $i | cut -d "/" -f 3)
TAXA=$(echo $i | cut -d "/" -f 3 | cut -d "." -f 1)

echo working on ${i}
# check that probes mapped correctly
echo "checking for map in bed file"
awk '$1 == $4 {print}' ../bed_files/${FASTA}.mapped.bed \
	> ../bed_files/${FASTA}.mapped.filt.bed

# make flank length as long as the file can be (generally) to allow for variation in sequence position
echo "flank"
bedtools flank -i ../bed_files/${FASTA}.mapped.filt.bed \
	-g ${i}.fai \
	-b 1200 \
	> ../bed_files/${FASTA}.flankloc.bed 

# extract flanking region.
echo "get fasta"
# extract the flanking regions
bedtools getfasta -nameOnly -fullHeader \
	-fi ${i} \
	-bed ../bed_files/${FASTA}.flankloc.bed \
	-fo ./${TAXA}.tmp.fasta

awk '/^>/ {
		if (prev == "") {
			prev = $0;
			printf("%s |%s\n", $0, substr($0, 2));
		} else if (prev == $0) {
			printf("--------------------"); # consider changing to "--------------------" so as not to break phyluces seq_cap_align ambigous bases filter
		} else {
			printf("\n%s |%s\n", $0, substr($0, 2));
			prev = $0;
		}
		next;
	}
	{
		printf("%s", $0);
	}
	END {
		printf("\n");
	}' ./${TAXA}.tmp.fasta | cut -d "_" -f 1,2,3 \
	> ../output//${TAXA}.flanks.fasta

echo -e "\n"
done

cd ../
rm -r tmp
