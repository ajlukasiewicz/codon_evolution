#!/usr/bin/env python3
import csv
import math

import numpy as np
import pandas as pd
from scipy import stats

#-------------------------------------------------------------------------------

# Imports table from pinetree output and calculates protein production rate
# at steady state
class EvalSteadyState:
    def __init__(self):
        self.time = []
        self.ribosome_counts = []
        self.protein_counts = []
        self.fitness = ""
        self.slopes = []
        self.x = []
        self.y = []

    def importTable(self,file,protein_name):
        '''open tsv file from pinetree, split each column into lists '''
        simulation = open(file, 'r')
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
        for i in range(len(feature)):
            if feature[i] == '__ribosome':
                self.ribosome_counts.append(proteins[i])
                self.time.append(time[i])
            if feature[i] == str(protein_name):
                self.protein_counts.append(proteins[i])
        return(self.time,self.ribosome_counts,self.protein_counts)

    def ribo_slope(self,times,ribosome_counts):
        for i in range(len(times)):
            if i <= len(times) - 3:
                self.slopes.append(float(ribosome_counts[i+1]
                                        - ribosome_counts[i])
                                        /(times[i+1]-times[i]))
        return self.slopes

    def steady_state(self,times,slopes,ribosome_counts,protein_counts):
        '''find where the slope of ribosome counts has leveled at zero,
            calculate the protein production rate after that point'''
        for i in range(len(times)):
            average_slope = np.average(self.slopes[i:i+6]) #average of 5 points in list
            if math.fabs(average_slope) <= 0.01 and i <= (len(times)-10): #finds slopes that were calculated closely around 0
                self.y = protein_counts[i:len(protein_counts)]
                self.x = times[i:len(protein_counts)]
                break
        return self.x, self.y

    def linreg(self,x,y):
        slope, intercept , r_value, p_value, std_err = stats.linregress(x,y)
        self.fitness = slope
        return self.fitness
