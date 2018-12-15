grep -e "Network" annotation.tsv | head -n1 > anno_all.tsv
grep -v -e "Network" annotation.tsv >> anno_all.tsv

sort anno_all.tsv | uniq > anno_dedup.tsv

grep -e "run" evaluation.tsv | head -n1 > eval_all.tsv
grep -v -e "run" evaluation.tsv >> eval_all.tsv

sort eval_all.tsv | uniq > eval_dedup.tsv
