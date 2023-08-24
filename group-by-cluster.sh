#!/bin/bash

source /local/anaconda3/bin/activate
conda activate characterization

GEN=$(echo $1 | cut -d . -f 1)
PROBE=$(echo $2 | cut -d . -f 1)
BP=$3

# get intergenic probes (we already have the intersect of genic probes)
echo -e "Find intersect of $PROBE and $GEN intergenic features \n"
bedtools intersect -a '1genefeatures/'${GEN}'.sorted.intergenic.gff' \
    -b '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bed' \
    -names ${PROBE} \
    -wb > '6group-by/'${PROBE}-${GEN}'.intergenic.intersect'

# get bed files of intergenic and genic loci
awk -F"\t" 'BEGIN{OFS="\t"} {print $10,$11,$12,$13,$14,$15}' '6group-by/'${PROBE}'-'${GEN}'.gene.intersect' > 6group-by/tmp.gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} {print $10,$11,$12,$13,$14,$15}' '6group-by/'${PROBE}'-'${GEN}'.intergenic.intersect' > 6group-by/tmp.intergenic.bed

# find genic and intergenic regions that intersect
bedtools intersect -a 6group-by/tmp.intergenic.bed \
	-b 6group-by/tmp.gene.bed | 
	awk '{print $4}' | \
	cut -d _ -f 1 | \
	sort -u \
	> 6group-by/genic-intergenic_loci.list

# get cluster from bedfile
bedtools cluster -i '6group-by/'${PROBE}'-'${GEN}'.mem.sorted.mapped.bed' \
	-d $BP > 6group-by/probes.cluster

# exclude intergenic-genic loci from cluster list
grep -vf 6group-by/genic-intergenic_loci.list 6group-by/probes.cluster > 6group-by/probes.filtered.cluster

# get bedtools cluster file grouped by the cluster
bedtools groupby -g 7 -c 4 -o collapse \
	-i 6group-by/probes.filtered.cluster \
	> 6group-by/tmp.probes.grpby_cluster

awk 'BEGIN{FS=OFS="\t"} {split($2, probes, ","); 
	for (i in probes) {split(probes[i], loci, "_"); 
	loci_array[loci[1]] = 1} delete probes; 
	for (j in loci_array)  {loci_list = loci_list ? loci_list "," j : j}; 
	delete loci_array; 
	print "cluster_"$1, loci_list; loci_list=""}' 6group-by/tmp.probes.grpby_cluster > '6group-by/'${PROBE}'-'${GEN}'.grpby_cluster'

