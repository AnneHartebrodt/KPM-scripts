#!/bin/#!/usr/bin/env bash

for file in $(ls $1)
do
java -jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/target/KPM-5-jar-with-dependencies.jar /home/anne/Documents/Master/MA/code/keypathwayminer-standalone/src/main/resources/batch_parameter_test.txt $1/$file /home/anne/Documents/Master/MA/Testing/log.txt
done
