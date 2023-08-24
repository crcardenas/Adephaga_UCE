---
title: Adephaga_UCE
author: CR Cardenas
date: 2023.08.217
---

# Adephaga_UCE README 

## UCE characerization and concatenation for phylogenomic analysis of Adephaga

The goal is to map probes to genomes in order to realize the total overalp of all probes used n the Adephaga data for __concatenation/merging multiple UCEs__ (as in [Van Dam et al 2021](https://academic.oup.com/sysbio/article/70/2/307/5880562)) on the same gene and __UCE characterization__. In general, folks will only have one probeset to use. But because I am integrating [anchored hybrid enrichment data](https://resjournals.onlinelibrary.wiley.com/doi/full/10.1111/syen.12508), [Adephaga UCE probes](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5260)[Adephaga UCE probes](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5260), and the original [Coleoptera UCE probes](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12754) I will need to create a merged dataset (nameed joined).

Two important assumptions being made about the genome and genes being used:
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

To follow this workflow, for now see the complete [markdown file](UCE_characterization_and_concatenation_workflow-v0.3.md).
