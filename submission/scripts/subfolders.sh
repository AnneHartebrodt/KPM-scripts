execute_R(){
  Rscript ~/Masterarbeit/R/reports/create_table_evaluation.R -f . -n /home/anne/Masterarbeit/finals/sample_networks/$net/ -g /home/anne/Masterarbeit/finals/sample_data/$net/ -o .
  cat annotation.tsv >> $1/annotation.tsv
  cat evaluation.tsv >> $1/evaluation.tsv
  cat fdr.tsv >> $1/fdr.tsv
}

net=$2 # biogrid
# $1=/home/anne/Masterarbeit/finals/result/
for m in greedy #random
do
  mkdir $m
  mv *_$m* $m
  cd $m
  for p in degreeaware nodeswap #edgerewire
  do
    mkdir $p
    mv *_$p* $p
    cd $p
    for t in sliding largest maximum
    do
      mkdir $t
      mv *_$t* $t
      cd $t
      for st in sum normDegSum meanLog # #
      do
      mkdir $st
      mv *_$st* $st
      cd $st
        for wt in 15 120 50 70
        do
          mkdir $wt
          mv *_$wt $wt
          cd $wt
          Rscript ~/Masterarbeit/R/reports/create_table_evaluation.R -f . -n /home/anne/Masterarbeit/final/sample_networks/$net/ -g /home/anne/Masterarbeit/final/sample_data/$net/ -o .
          cat annotation.tsv >> $1/annotation.tsv
          cat evaluation.tsv >> $1/evaluation.tsv
          cd ..
        done
        cd ..
      done
      cd ..
    done
    cd ..
done
cd ..
done
