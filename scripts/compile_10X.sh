#!/bin/bash

display_usage() { 
        echo -e "\nUsage:\n$0 [-i input_file] [-o output_file]\n" >&2
}

while getopts "hi:o:s:r:" opt; do
        case "$opt" in
		h)	display_usage
			exit 1;;
                i)	input_file="$OPTARG" ;; # (str);;
                o)	output_file="$OPTARG" ;; # e.g. batch1
		s)	sample_id="$OPTARG" ;;
		r)	ref_id="$OPTARG" ;;
        esac
done


consensusString=$(sed -n '2p' ${input_file})
consensusBases=$(sed 's/N//g' <<< ${consensusString})
consensusBasePercent=$(echo "scale=4;${#consensusBases} / ${#consensusString}" | bc)
printf "%s\t%d\t%d\t%.4f\n" ${sample_id}_${ref_id} ${#consensusString} ${#consensusBases} ${consensusBasePercent} >> ${output_file}
#printf "%s\t%d\t%d\t%.4f\n" !{sample_id}_!{ref_id} ${consensusString} ${consensusBases} ${consensusBasePercent} >> compiled_10x_coverage.tsv

