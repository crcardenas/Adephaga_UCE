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
	
	echo "Pairwise comparison between $species_1 and $species_2 completed. Output saved to $newdirectory/$output_file"

	done

done