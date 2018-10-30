#!/bin/bash
for fn in $(cut -f11 $1 | cut -f1 -d";")
do
echo $fn
wget $fn
done 

