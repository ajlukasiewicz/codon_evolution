#!/usr/bin/env python

import csv
import numpy as np
import pandas as pd
import math
import scipy
from scipy import stats
#-------------------------------------------------------------------------------

simulation = open("function_test_counts.tsv", 'r')
proteins = []
time = []
feature = []
lines = simulation.readlines()
linecounter = 0
for line in lines:
    if linecounter == 0:
        print("header found")
    else:
        data = line.strip().split('\t')
        proteins.append(float(data[2]))
        time.append(float(data[0]))
        feature.append(str(data[1]))
    linecounter += 1
simulation.close()

#filter table for protein feature
ribosome_counts = []
protein_counts = []
times = []
for i in range(len(feature)):
    if feature[i] =='__ribosome':
        ribosome_counts.append(proteins[i])
        times.append(time[i])
    if feature[i] == 'proteinP':
        protein_counts.append(proteins[i])

list_len = len(times)
#print(list_len)

#Calculate slope based on feature counts over time
slopes = []
for i in range(len(times)):
    if i <= list_len - 3:
        slopes.append(float(ribosome_counts[i+1] - ribosome_counts[i])/(times[i+1]-times[i]))

#write output table
production_rates = ['Time,Ribosomes,Proteins,ProductionRate\n']

#find steady state production rates
for i in range(len(times)):
    average_slope = np.average(slopes[i:i+6]) #average of 5 points in list
    if math.fabs(average_slope) <= 0.01 and i <= (len(times)-10): #finds slopes that were calculated closely around 0
        print("time: %.2f " % times[i])
        print("free ribosomes: %.2f " % ribosome_counts[i])
        print("amount of proteinP produced: %.2f" % protein_counts[i])
        production_rate = (protein_counts[i+6] - protein_counts[i])/(times[i+6]-times[i])
        print("steady state production rate calculated: %.2f" % production_rate)
        production_rates.append(str(times[i])+','+ str(ribosome_counts[i]) + ',' + str(protein_counts[i]) +','+ str(production_rate)+ '\n')

with open('function_test_production_rates.csv', "w") as f:
    f.writelines(production_rates)
f.close()
