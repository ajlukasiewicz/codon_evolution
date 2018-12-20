#!/usr/bin/env python

import argparse
import csv
import math
import numpy as np
import pandas as pd
from scipy import stats
#-------------------------------------------------------------------------------

def slope(times,slopes,ribosome_counts):
    for i in range(len(times)):
        if i <= len(times) - 3:
            slopes.append(float(ribosome_counts[i+1] - ribosome_counts[i])/(times[i+1]-times[i]))
    return slopes

def steady_state(times,slopes,ribosome_counts,protein_counts,x,y):
    for i in range(len(times)):
        average_slope = np.average(slopes[i:i+6]) #average of 5 points in list
        if math.fabs(average_slope) <= 0.01 and i <= (len(times)-10): #finds slopes that were calculated closely around 0
            print("time: %.2f " % times[i])
            print("free ribosomes: %.2f " % ribosome_counts[i])
            print("amount of proteinP produced: %.2f" % protein_counts[i])
            y = protein_counts[i:len(protein_counts)]
            x = times[i:len(protein_counts)]
            break
    slope, intercept , r_value, p_value, std_err = stats.linregress(x,y)
    print("steady state production rate calculated: %.2f" % slope + " proteins/s")
    return slope

def main():
    parser = argparse.ArgumentParser(description='simulation counts table')
    parser.add_argument(
        '-i',
        action='store',
        dest='i',
        required=True,
        type=str,
        help="input .tsv file from pinetree",
        )
    #open tsv file from pinetree, split each column into lists
    options = parser.parse_args()
    simulation = open(options.i, 'r')
    proteins = []
    time = []
    feature = []
    lines = simulation.readlines()
    linecounter = 0
    for line in lines:
        if linecounter == 0:
            continue
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

    #calculate slope of free ribosome counts
    slopes = []
    slope(times,slopes,ribosome_counts)

    #determine time and rate of protein production at steady state
    x = []
    y = []
    steady_state(times,slopes,ribosome_counts,protein_counts,x,y)

if __name__ == "__main__":
    main()
