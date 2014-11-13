#!/bin/bash


for j in 1 2 3 4 5
do
    for i in 0.1 0.05 0.15 0.2 0.25
    do
        EvolveAGene4 -f hxb2_env -t ran -b $i -ss both -n 32 -a 0.2
        ramsaics -i hxb2_env_Unaligned.FASTA -o out_noah.fasta --clean --arc-align > /dev/null 2>&1
        qscore -test out_noah.fasta -ref hxb2_env_True_alignment.FASTA >> noah_scores.txt
        ramsaics -i hxb2_env_Unaligned.FASTA -o out_noah_no_fb.fasta --clean > /dev/null 2>&1
        qscore -test out_noah_no_fb.fasta -ref hxb2_env_True_alignment.FASTA >> noah_no_fb_scores.txt
        ramsaics -i hxb2_env_Unaligned.FASTA -o out_noah_basic.fasta --clean --no-arc-guide > /dev/null 2>&1
        qscore -test out_noah_basic.fasta -ref hxb2_env_True_alignment.FASTA >> noah_basic_scores.txt
        ramsaics -i hxb2_env_Unaligned.FASTA -o out_noah_no_arc.fasta --clean --arc-align --no-arc-guide > /dev/null 2>&1
        qscore -test out_noah_no_arc.fasta -ref hxb2_env_True_alignment.FASTA >> noah_no_arc_scores.txt
        muscle -in hxb2_env_Unaligned.FASTA -out out_muscle.fasta -maxiters 2 > /dev/null 2>&1
        qscore -test out_muscle.fasta -ref hxb2_env_True_alignment.FASTA >> muscle_scores.txt
        clustalo -i hxb2_env_Unaligned.FASTA -o out_clustalo.fasta --force > /dev/null 2>&1
        qscore -test out_clustalo.fasta -ref hxb2_env_True_alignment.FASTA >> clustalo_scores.txt
        clustalw -infile=hxb2_env_Unaligned.FASTA -outfile=out_clustalw_temp.fasta > /dev/null 2>&1
        seqret -sequence=out_clustalw_temp.fasta -outseq=out_clustalw.fasta -osformat2 fasta > /dev/null 2>&1
        qscore -test out_clustalw.fasta -ref hxb2_env_True_alignment.FASTA >> clustalw_scores.txt
        mafft --auto hxb2_env_Unaligned.FASTA > out_mafft.fasta
        qscore -test out_mafft.fasta -ref hxb2_env_True_alignment.FASTA -ignoretestcase -ignorerefcase >> mafft_scores.txt
        prank +F -d=hxb2_env_Unaligned.FASTA -o=out_prank_f.fasta -iterate=1 -codon > /dev/null 2>&1
        qscore -test out_prank_f.fasta.best.fas -ref hxb2_env_True_alignment.FASTA >> prank_f_scores.txt
        prank -d=hxb2_env_Unaligned.FASTA -o=out_prank.fasta -iterate=1 -codon > /dev/null 2>&1
        qscore -test out_prank.fasta.best.fas -ref hxb2_env_True_alignment.FASTA >> prank_scores.txt
    done
done

