#!/usr/bin/env bash
find . -name 'annotation.tsv' >> anno.tsv
while read line; do $(cat "$line" >> file.tsv); done < anno.tsv
grep -e "Network" file.tsv | head -n1 > anno_all.tsv
grep -v -e "Network" file.tsv >> anno_all.tsv

find . -name evaluation.tsv >> eval.tsv
while read line; do $(cat "$line" >> file2.tsv); done < eval.tsv
grep -e "run" file2.tsv | head -n1 > eval_all.tsv
grep -v -e "run" file2.tsv >> eval_all.tsv

find . -name fdr.tsv >> f.tsv
while read line; do $(cat "$line" >> file3.tsv); done < f.tsv
grep -e "V1" file3.tsv | head -n1 > fdr_all.tsv
grep -v -e "V2" file3.tsv >> fdr_all.tsv
