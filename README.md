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


## Data

This pipeline will mainly use the [*Pterostichus madidus*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9315234/) genome for testing and finalizing the workflow, but others will be used. Genomes with annotations (general feature files, GFF) can be found at the [Darwin Tree of Life ensembl](https://projects.ensembl.org/darwin-tree-of-life/). 

Ensure you have a GFF3 format. For details on GFF format, see [here](http://gmod.org/wiki/GFF3).

### __*!!! genomes used are likely to change by final analysis, likely just the ones with GFF files Neb.brevicolis, Oph.ardosiacus, and Pte.madidus*__

__Available genomes and their names names__

|code|name|ID|
|---|---|---|
|agrPla1|*Agrilus planipennis*|GCA_000699045.1|
|\* ampIns|*Amphizoa insolens*|DNA3784|
|anoGla1|*Anoplophora glabripennis*|GCA_000390285.1|
|\* bemHap1|*Bembidion haplogonum*|DNA2544|
|\* chlSer1|*Chlaenius sericeus*|DNA4821|
|denPon1|*Dendroctonus ponderosae*|GCA_000355655.1|
|lepDec|*Leptinotarsa decemlineata*|GCA_000500325.1|
|lioTuu1|*Lionepha*|DNA3782|
|menMol1|*Mengenilla moldrzyki*|GCA_000281935.1|
|! nebBrevi1|*Nebria brevicolis*|GCA_944738965.1|
|! nebSali1|*Nebria salina*|GCA_944039245.1|
|\* omoHam1|*Omoglymmius hamatus*|DNA3783|
|ontTau1|*Onthophagus taurus*|GCA_000648695.1|
|! ophArdo1|*Ophonus ardosiacus*|GCA_943142095.1|
|! pteMad2|*Pterostichus madidus*|GCA_911728475.2|
|\* pteMel1|*Pterostichus melenarius*|DNA3787|
|\* traGib1|*Trachypachus insolens*|DNA3786|
|triCas2|*Tribolium castaneum*|GCA_000002335.2|

(* indicates Adephaga reference genome; ! indicates a genome feature file)

__probe file names:__
* vasilikopoulos_etal_ahe_probes.fasta
* Adephaga_2.9Kv1_UCE-Probes.fasta
* Coleoptera-UCE-1.1K-v1.fasta

### genomes are removed from the github repository due to their size

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

This pipeline should be able to run on a personal laptop (windows with linux subsystem and linux, uncertain about mac) with sufficient storage, memory and CPU available. Alternatively, you can .... [add json information]
