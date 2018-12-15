#!/bin/bash
#for fn in fastq/SRR17472{08..11};
#do
#samp=`realpath ${fn}`
#echo "Processing sample ${samp}"
#salmon quant -i ./genomes/human/human_index -l A -1 ${samp}_1.fastq.gz -2 ${samp}_2.fastq.gz \
#         -p 3 -o quants/${fn}_quant --gcBias 
#done
#
#mkdir quants
for fn in SRR1747{184..211};
do
samp=`realpath ${fn}`
echo "Processing sample ${samp}"
salmon quant -i ../genomes/human/human_index -l A -1 ${samp}_1.fastq.gz -2 ${samp}_2.fastq.gz \
         -p 3 -o quants/${fn}_quant --gcBias
done 