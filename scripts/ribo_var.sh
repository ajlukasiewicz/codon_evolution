#!/bin/bash

for ((i=5;i<=20; i++))
do
	python3 pooled_evolution.py -g 1000 -r"$i" -o ribosome_num_"$i" 
done
