#!/usr/bin/env python

import csv
import numpy as np
import pandas as pd
import math

#-------------------------------------------------------------------------------

def slope(times,slopes,ribosome_counts):
    for i in range(len(times)):
        if i <= len(times) - 3:
            slopes.append(float(ribosome_counts[i+1] - ribosome_counts[i])/(times[i+1]-times[i]))

def steady_state(times,slopes,ribosome_counts,protein_counts):
    production_rates = ['Time,Ribosomes,Proteins,ProductionRate\n']
    for i in range(len(times)):
        average_slope = np.average(slopes[i:i+6]) #average of 5 points in list
        if math.fabs(average_slope) <= 0.01 and i <= (len(times)-10): #finds slopes that were calculated closely around 0
            print("time: %.2f " % times[i])
            print("free ribosomes: %.2f " % ribosome_counts[i])
            print("amount of proteinP produced: %.2f" % protein_counts[i])
            production_rate = (protein_counts[i+6] - protein_counts[i])/(times[i+6]-times[i])
            print("steady state production rate calculated: %.2f" % production_rate)
            production_rates.append(str(times[i])+','+ str(ribosome_counts[i]) + ',' + str(protein_counts[i]) +','+ str(production_rate)+ '\n')
