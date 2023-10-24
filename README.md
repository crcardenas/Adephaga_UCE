---
title: Adephaga_UCE
author: CR Cardenas
reviewed by: Dr. Jeremy Gauthier
date: 2023.09.26
---

# Adephaga_UCE README 

## Tasks

1. double check synteny workflow
2. Add slice flanking workflow
3. add R markdown/script files
4. check grammar in markdown files

## UCE characerization and concatenation for phylogenomic analysis of Adephaga

The goal is to map probes to genomes in order to realize the total overalp of all probes used n the Adephaga data for __concatenation/merging multiple UCEs__ (as in [Van Dam et al 2021](https://academic.oup.com/sysbio/article/70/2/307/5880562)) on the same gene and __UCE characterization__. In general, folks will only have one probeset to use. But because I am integrating [anchored hybrid enrichment data](https://resjournals.onlinelibrary.wiley.com/doi/full/10.1111/syen.12508), [Adephaga UCE probes](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5260)[Adephaga UCE probes](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5260), and the original [Coleoptera UCE probes](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12754) I will need to create a merged dataset (nameed joined).

Two **important assumptions** being made about the genome and genes being used:
1. Intergenic sequences close to genes are not being considered as promotors. This is mainly due to this information being unknown.
2. There is no alternative splicing or overlapping genes (ex: [overlaping genes](https://doi.org/10.1186/s12864-021-08039-6)).


## Goals of this workflow:
1. identify what probes are genetic or intergenic 
2. identify the overlap/intersect of probesets used between datasets
3. create a new probeset that integrates all probes for use in [Phyluce](https://phyluce.readthedocs.io/en/latest/)
4. create a list of probes that should be concatenated in a partition file for phylogenetic analysis
    1. create a script to integrate it based on an existing partition file (e.g., output of Phyluce)

# Software and packages used

### __*!!! create a json file with environment for download and easy duplication of the environent*__

* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
* [BWA](https://bio-bwa.sourceforge.net/)
* [natsort](https://anaconda.org/anaconda/natsort)
* [sam tools](http://www.htslib.org/)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [R](https://www.r-project.org/)

Recomended installation procedure:
```
conda create --name characterization
conda activate characterization
conda install -c bioconda bedtools blast bwa samtools seqkit 
conda install -c anaconda natsort
```

This pipeline should be able to run on a personal laptop (windows with linux subsystem and linux, uncertain about mac) with sufficient storage, memory and CPU available.

All scripts are found in the [scripts directory](/scripts/)

To follow this workflow, for now see the complete [markdown files](/markdown_files/workflow_v4_26092023.md); or follow the step by step workflow:
1. [Data processing](/markdown_files/data_processing.md)
2. [Map probes to gene features and intersect](/markdown_files/map_probes_and_intersect.md)
3. [Integrate probesets](/markdown_files/integrate_pobesets.md)
4. [Concatenate probes by gene](/markdown_files/concatenate_loci.md)

Lastly, to see the synteny comparison of genomes find the [synteny workflow here](/synteny_workflow/synteny_workflow.md)

## Additional Adephaga UCE scripts

An additional downstream result of this workflow is testing the effect of flanking+core-probe-regions, flanking, and core-probe-regions of sequence data in phylogenetic analysis inside phyluce. Initial analyses showed incredibly gappy alignments, likely due to the depth of evolutionary relationships and integration of a diverse set of genomic data collection methods. 

The following markdown file: contains the scripts necessary to "slice" flanking regions from UCE data.

1. [Slice flanking from core probe region](https://github.com/crcardenas/Adephaga_UCE/blob/main/markdown_files/Slice_flanking_from_core.md)
