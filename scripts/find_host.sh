#!/bin/bash

SAMPLE=$1


var=$(find /mnt/phage_01/MIDRP/references/ -maxdepth 1 -name "*.fasta" | xargs -I{} basename {} | awk -F_ '{print $1FS$2}')

for x in ${var}
do
	#echo $x
	if [[ "${SAMPLE}" =~ "${x}" ]]
	then
		reffasta=$(find /mnt/phage_01/MIDRP/references/ -maxdepth 1 -name "*${x}*")
		echo "$reffasta"
		break
	fi
done
