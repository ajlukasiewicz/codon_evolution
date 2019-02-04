#!/bin/bash

for ((i=1; i<=100; i++))
do
    python evolve_var_rate.py -g1 -r"$i" -m 100 &  

done