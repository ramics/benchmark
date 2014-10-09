#!/bin/bash

declare -a tools=("noah_full" "prank" "prank_f")
declare -a bali=("11" "12" "20" "30" "40" "50")

for i in "${bali[@]}"
do 
        for tool in "${tools[@]}"
	do
		python run.py -i bb3_release/RV${i}/ -o out_${i} -c scripts/tools/${tool}.sh 
	done
done
