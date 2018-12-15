#!/bin/bash

outdir="/home/h/hartebrodt"
mkdir "$outdir/reports"
mkdir "$outdir/quants"

for fn in $(cut -f11 $1)
do

out=$(basename $fn | cut -d. -f1 | cut -d_ -f1)
echo $out

file1=$(echo $fn | cut -f1 -d";")
file2=$(echo $fn | cut -f2 -d";")
#echo $file1
#echo $file2

wget $file1
wget $file2
samp1=$(basename $file1)
samp2=$(basename $file2)
echo "Processing sample ${out}"

samp1out=out.$samp1
samp2out=out.$samp2



echo "/home/h/hartebrodt/fastp -i $samp1 -I $samp2 -o $samp1out -O out.$samp2out -q 20 --disable_length_filtering --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -h=reports/$out.html"

rm $samp1
rm $samp2


echo "/home/h/hartebrodt/salmon-0.11.3-linux_x86_64/bin/salmon quant -i /big/h/hartebrodt/rnaseq/genome/human_index -l A -1 $samp1out -2 $samp2out \
         -p 3 -o $outdir/quants/${out}_quant --gcBias"


done 

