#!/bin/bash

for i in {1..5}
do
    python generate.py NewNamesNucSeqs.fasta > /dev/null 2>&1
    time ramsaics -i NewNamesNucSeqs.fasta.in.fasta -o out_ramsaics.fasta --no-arc-guide --arc-align --roche > /dev/null
    qscore -ref NewNamesNucSeqs.fasta.out.fasta -test out_ramsaics.fasta -ignorerefcase -ignoretestcase >> noah_scores.txt
    time java -jar ~/Downloads/macse_v1.01b.jar -prog alignSequences -fs 10 -stop 15 -seq ~/Documents/Thesis/Code/benchmarking/benchmark/frameshifts/NewNamesNucSeqs.fasta.in.fasta  -out_NT ~/Documents/Thesis/Code/benchmarking/benchmark/frameshifts/out_macse.fasta > /dev/null 2>&1
    sed 's/!/-/g' out_macse.fasta > out_macse_final.fasta
    qscore -ref NewNamesNucSeqs.fasta.out.fasta -test out_macse_final.fasta -ignorerefcase -ignoretestcase >> macse_scores.txt
    time mafft --localpair --maxiterate 2 NewNamesNucSeqs.fasta.in.fasta > out_mafft.fasta 2>/dev/null
    qscore -ref NewNamesNucSeqs.fasta.out.fasta -test out_mafft.fasta -ignorerefcase -ignoretestcase >> mafft_scores.txt
    time clustalo -i NewNamesNucSeqs.fasta.in.fasta -o out_clustalo.fasta --force
    qscore -ref NewNamesNucSeqs.fasta.out.fasta -test out_clustalo.fasta -ignorerefcase -ignoretestcase >> clustalo_scores.txt

done
