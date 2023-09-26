 
#!/bin/bash

species=("nebBrev1" "nebSali1" "ophArdo1" "pteMadi2")

for ((i=0; i<${#species[@]}; i++)); do
    for ((j=i+1; j<${#species[@]}; j++)); do

		species_1="${species[i]}"
    	species_2="${species[j]}"

		# preprocessing steps

		Rscript --vanilla script_synteny_v2.R "${species_1}" "${species_2}"

	done
done
