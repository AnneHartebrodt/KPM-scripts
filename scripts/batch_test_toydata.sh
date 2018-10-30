#!/bin/#!/usr/bin/env bash

## declare an array variable
declare -a arr=("mean" "median" "sum")

## now loop through the above array
for method in "${arr[@]}"
do
for file in $(ls | grep "samplewise")
do
echo $file
file2="${file/samplewise/general}"
outdir="${file/toy_samplewise_/""}"
outdir="${outdir/\.tsv/""}"
outdir=$outdir\/$method
echo $outdir

#-numProc=1 -matrix1=$file -datasetsFile=/home/anne/Documents/Master/MA/Testing/datasets.txt -summaryFile=harkan1.txt -combineOp=OR -pathwaysStatsFile=harkan2.txt -resultsDir=$outdir -graphFile=/home/anne/Documents/Master/MA/data/networks/Biog.sif -geneStatsFile=harkan3.txt -K=5 -L1=2 -maxsolutions=5 -perturbation_technique=edgerewire -algo=FDR -strategy=FDR -Umove_bens -comparator=LET -significance_level=0.05 -L1_pvalues=$file2 -use_double -mfHeader -validation_file=/home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/COAD-VAL-ENTREZ.txt -L1_pvaluecutoff=0.05 -randomized_graph_file=/home/anne/Documents/Master/MA/data/networks/Biog_randomized.sif -aggregation_method=$method

Rscript ~/Masterarbeit/R/reports/create_report_greedy.R -f $outdir -g ~/Masterarbeit/data/toydata/expression_type.tsv

done
done
