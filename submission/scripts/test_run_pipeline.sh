#!/bin/bash
#make Test networks
#bash ~/Masterarbeit/scripts/batch_exec_KPM.sh ~/Masterarbeit/test_runs/generate_samples/

#generate Testdata with different number of networks hidden in data
for i in {1..3}
do
Rscript ../R/real_toydata.R -f ../differential_out/dataHD/DESEQ/results.tsv -n "String900" -m Homo_Sapiens_String_NCBI900.tsv -b . -k ../data/networks/StringDB/ -r 5 -d $i
Rscript ../R/real_toydata.R -f ../differential_out/dataHD/DESEQ/results.tsv -n "String700" -m Homo_Sapiens_String_NCBI700.tsv -b . -k ../data/networks/StringDB/ -r 5 -d $i
Rscript ../R/real_toydata.R -f ../differential_out/dataHD/DESEQ/results.tsv -n "biogrid" -m Homo_sapiens_biogrid.tsv -b . -k ../data/networks/biogrid/ -r 5 -d $i
done






