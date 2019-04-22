#!/bin/bash

for ((i=5;i<=6; i++))
do
	mkdir -p ../data/ribosome_num_"$i"
	python3 pooled_evolution.py -g 1 -s 0.5 -f 1.0 -r"$i"
	mmv -m ../data/*.tsv ribosome_num_"$i"
	mmv -m ../data/*.csv ribosome_num_"$i"
done
