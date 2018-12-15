#!/bin/bash
#Rscript /home/anne/Masterarbeit/R/differential/differential_gene_expression_HD.R -q /home/anne/Documents/Master/MA/data/dataHD/quants/ -m /home/anne/Documents/Master/MA/data/identifier_mappings/mart_export.txt -r "/home/anne/Documents/Master/MA/data/dataHD/SraRunTable(2).txt" -o /home/anne/Documents/Master/MA/differential_out/dataHD/DESEQ_final/ -t /home/anne/Documents/Master/MA/real_data/data/ -n /home/anne/Documents/Master/MA/data/networks/StringDB/Homo_Sapiens_String_NCBI900.tsv -p hsa05016 -a /home/anne/Documents/Master/MA/data/identifier_mappings/mart_export_names.txt

#Rscript /home/anne/Masterarbeit/R/bionet/bionet_test_1.R -o Masterarbeit/bionet_out/dataHD/StringDB/ -f /home/anne/Documents/Master/MA/differential_out/dataHD/DESEQ_final/results.tsv -n Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_ENS900.tsv

#Rscript /home/anne/Masterarbeit/R/bionet/bionet_results_eval.R -f Masterarbeit/bionet_out/dataHD/StringDB/module.sif -n Masterarbeit/data/networks/StringDB/Homo_Sapiens_String_ENS900.tsv -m Masterarbeit/data/identifier_mappings/mart_export.txt -p  hsa05016 -o Masterarbeit/bionet_out/dataHD/StringDB/

time java -jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/target/KPM-5-jar-with-dependencies.jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/batch_parameter_test.txt ~/Masterarbeit/real_data/call/call120_001.txt /home/anne/Documents/Master/MA/Testing/log.txt

time java -jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/target/KPM-5-jar-with-dependencies.jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/batch_parameter_test.txt ~/Masterarbeit/real_data/call/call_15_001.txt /home/anne/Documents/Master/MA/Testing/log.txt

time java -jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/target/KPM-5-jar-with-dependencies.jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/batch_parameter_test.txt ~/Masterarbeit/real_data/call/call_15_005.txt /home/anne/Documents/Master/MA/Testing/log.txt

time java -jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/target/KPM-5-jar-with-dependencies.jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/batch_parameter_test.txt ~/Masterarbeit/real_data/call/call_120_005.txt /home/anne/Documents/Master/MA/Testing/log.txt
