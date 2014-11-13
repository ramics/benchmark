#!/bin/bash

declare -a tools=("noah_full" "noah_basic" "noah_no_arc" "noah_no_fb")

for tool in "${tools[@]}"
do
	echo ${tool}
	time python run.py -r DNA_alignments -i random_selection -o out -c scripts/tools/${tool}.sh > /dev/null 2>&1 
done
