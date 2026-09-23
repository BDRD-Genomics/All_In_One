#!/bin/bash

all_parsed_out=$1

#using tail so we don't collect header in the awk command
tail -n +2 ${all_parsed_out} | awk -F"," '{print $2}' | sort -u > all_node_ids.txt
awk -F"," '$18 ~ /Viruses/ {print $2}' ${all_parsed_out} | sort -u > viral_node_ids.txt
comm -13 viral_node_ids.txt all_node_ids.txt > non-viral_node_ids.txt
#diff all_node_ids.txt viral_node_ids.txt | grep -E '^<' | cut -c3- > non-viral_node_ids.txt

#grab file sizes
declare -i non_viral_size
declare -i viral_size
non_viral_size=$(wc -c non-viral_node_ids.txt | cut -d' ' -f1 )
viral_size=$(wc -c viral_node_ids.txt | cut -d' ' -f1 )

echo "File size of viral_node_ids.txt is $viral_size"
echo "File size of non-viral_node_ids.txt is $non_viral_size"

#write to output
head -1 ${all_parsed_out} > putative_viral_blast.csv
if [[ $non_viral_size -gt $viral_size ]];then
	echo "using viral_node_ids.txt to filter putative_viral_blast.csv"
	grep -Ff viral_node_ids.txt ${all_parsed_out} >> putative_viral_blast.csv
else
	echo "using non-viral_node_ids.txt to filter putative_viral_blast.csv"
	grep -Fvf non-viral_node_ids.txt ${all_parsed_out} >> putative_viral_blast.csv
fi
