# Slice flanking from core probe regions

This pipeline is intented to investigate how gappy flanking regions capture methods impact phylogenetic reconstruction. However, inorder to do so those flanking regions need to be identified and extracted from loci. This pipeline is a part of the [Adephaga_UCE](https://github.com/crcardenas/Adephaga_UCE) workflow and depends on the phyluce pipeline. This is soley for simplicity and consistency with the Adephaga_UCE project.

### Overview

In preliminary untrimmed alignments it was found that many genomic loci captured with AHE and UCE alignmetns were extremely gappy and unexpectedly long. Some loci were upwards of 20k bp vs the expected < 1k bp seen in UCE datasets. To help reduce this effect, [`phyluce_probe_slice_sequence_from_genomes`](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-3.html#extracting-fasta-sequence-matching-uce-loci-from-genome-sequences) flanking parameter was shortend to 400 bp for all non-UCE data rather than the 800 for transcriptomes and 1600 for genomic data, as described by [Bossert et al, 2019](https://www.sciencedirect.com/science/article/pii/S1055790318304317#f0005). This somewhat helped reduced gappy-ness but loci still exceded 2-3kb lengths. In addition, it was found that some assembled UCE loci were exceeding this length as well. To help resolve this issue, all data types were trimmed to 400 bp. While the gapppy data still remained in all cases previosly desribed there often remained a core region that was well aligned, the targeted probe regions. Moreover, when comparing data trimmed with [gblocks](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#aligning-uce-loci) with untrimmed data as advised by the phyluce pipeline loci often were liberally trimmed; for example 1800 bp to 300 bp. Potentially removing parsimony informative sites. An additional trimal analyses was performed, in another step of the Adephaga_UCE workflow, but it is out of the scope of this document. While there are a few reasons for the gaps and how the loss of sites from GBLOCKS might effect phylogenetic analyses, discussion of this is further outside the scope of this workflow and should be expected to be discussed in the forthcoming article.

The following workflow uses the phyluce pipeline [Harvesting UCE loci from genomes](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-3.html) with the flanking modified to 400bp (core+flanking) and 0 bp (core) up to the [Finding UCE ](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#finding-uce-loci) step in the work flow. These two sets are used to ID the flanking regions and slice them from the loci. __Keep in mind__, this workflow deviates from typical UCE workflows with phyluce, because it uses an integrated probe set (AHE and two UCE probe sets).

In order to investigate the "gap issue," three data sets were generated:

1. core+flanking: probe and flanking region loci, e.g. those generated with phyluce, where all loci have up to 400 bp flanking regions
2. core: core probe region, what is targeted by capture probes
3. flanking: where the core probe region is removed from the probe leaving flanking regions


## Workflow 

__software requirements__: bedtools, mafft, phyluce, seqkit, bwa ...

phyluce_assembly_explode_get_fastas_file was run with the additional `--taxon` parameter to produce the `exploded-fastas-*` directories. 
Ensure that the slice_flanking structure is constructed with the `bed_files`, `loci`, `noflank`, and `output` directories.

```
Adephaga_2023
├──
...
├── contigs_joined
├── contigs_joined_noflank
├── uce-search-results-joined
└── uce-search-results-joined_noflank_max200bp
	taxon-sets
	├── joined
	│   ├── exploded-fastas-loci
	│   ├── exploded-fastas-taxon
	│   └── log
	├── joined_flanking
	│   ├── exploded-fastas-loci
	│   ├── exploded-fastas-taxon
	│   ├── log
	│   └── slice_flanking
	│       ├── bed_files
	│       ├── loci
	│       ├── noflank
	│       └── output
	└── joined_noflank_max200bp
		├── exploded-fastas-loci
		├── exploded-fastas-taxon
		└── log
```

Loci are filtered if longer than 200 bp in the no_flank_max bp dataset.
```
mkdir contigs_joined_noflank_max200bp/
for i in contigs_joined_noflank/*.fasta; do NEW=$(echo $i | cut -d "/" -f 2 | cut -d "." -f 1); seqkit seq -M 200 ${i} > contigs_joined_noflank_max200bp/${NEW}.fasta; done
```

First start by creating symlink of flanking+core `joined_subset/exploded-fasta-taxon/*.unaligned.fasta` into the `joined-subset_noflanking/slice_flanking/loci`; the joined_subset contains both the core and flanking regions of interest. Repeat this process for the core probe regions as well.

```
# be literal with pathways, otherwise symlinks break in my experience.
for i in $(ls joined-subset/exploded-fastas-taxon); do
ln -s /data/work/Calosoma_phylo/phylogeny/2023_adephaga/taxon-sets/joined-subset/exploded-fastas-taxon/${i} \
/data/work/Calosoma_phylo/phylogeny/2023_adephaga/taxon-sets/joined-subset_flanking/slice_flanking/loci/${i}; 
done

for i in $(ls joined-subset_noflank-max200bp/exploded-fastas-taxon/); do 
ln -s /data/work/Calosoma_phylo/phylogeny/2023_adephaga/taxon-sets/joined-subset_noflank-max200bp/exploded-fastas-taxon/${i} \
/data/work/Calosoma_phylo/phylogeny/2023_adephaga/taxon-sets/joined-subset_flanking/slice_flanking/noflank/${i}; 
done
```

Next, we use the `slice_flanking.sh` script in the `joined-subset_flanking/slice_flanking/` directory. This  script first indexes all the core+flanking unaligned loci. Then maps the core loci against the core+flank loci using [BWA mem](https://bio-bwa.sourceforge.net/bwa.shtml). Next the bed file is checked to ensure that the mapped loci are the same and up to 1200 bp of flanking length is used (flanking + core: 400 bp + core bp + 400 bp loci are generally no greater than 1200 bp). Lastly, bedtools extracts the loci and an awk script pastes the flanking and 20 ambigous empty bases between each flanking region.
```
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
```

Lastly, create a file of all fastas for phyluce pipeline
```
cat output/*.fasta > ../all-taxa-incomplete.fasta
```

The phyluce workflow can be continued with the `phyluce_assembly_explode_get_fastas_file` step in the `joined-subset_flanking` directory. 

```
source /local/anaconda3/bin/activate
conda activate phyluce-1.7.1

phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas-taxon \
    --by-taxon

phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas-uce

```

---

Note: If using `phyluce_align_seqcap_align` for alignments and you choose to change line 70 from `-` to `N` then add --ambigous command in pipeline.
