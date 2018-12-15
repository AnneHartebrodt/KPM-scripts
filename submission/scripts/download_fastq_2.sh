#!/bin/bash
#move to target directory $2
#$1 = European Read Archive Run selector file
cd $2
for fn in $(cut -f11 $1 | cut -f2 -d";")
do
echo $fn
wget $fn
done 
# Download the fasta files of the RNA Seq experiments
