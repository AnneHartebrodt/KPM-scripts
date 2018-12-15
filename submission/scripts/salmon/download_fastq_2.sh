#!/bin/bash
for fn in $(cut -f11 $1 | cut -f2 -d";")
do
echo $fn
wget $fn
done 