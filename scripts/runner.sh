#!/bin/#!/usr/bin/env bash

#bash ../scripts/batch_exec_KPM.sh call/biogrid_sum/
#bash ../scripts/batch_exec_KPM.sh call/biogrid_normSum/
#bash ../scripts/batch_exec_KPM.sh call/biogrid_normDegSum/
#bash ../scripts/batch_exec_KPM.sh call/StringDB_sum/
#bash ../scripts/batch_exec_KPM.sh call/StringDB_normSum/


for file in $(ls $1)
do
bash ~/Masterarbeit/scripts/batch_exec_KPM.sh $1/$file
done