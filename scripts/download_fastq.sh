#!/bin/bash
#move to target directory
#$1 = European Read Archive Run selector file
cd $2
for fn in $(cut -f11 $1 | cut -f1 -d";")
do
echo $fn
wget $fn
done 

