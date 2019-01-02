#!/usr/bin/env python
'''
Define codon translation rates for use in pinetree simulation
Author: Alexandra Lukasiewicz
'''
import numpy as np

#-------------------------------------------------------------------------------
class Transcript:
#include transcript features needed for pt here? (example: if adding an additional gene, or changing the binding site)
    rates = [0.5,1.0]

    def __init__(self,length,generation):
        self.length = length
        self.weights = []
        self.generation = generation

    def random_codons(self):
        for i in range(0,self.length):
            r = np.random.randint(0,2)
            self.weights.append(self.rates[r])
        return self.weights

    def defined_codons(self,rate):
        for i in range(0,self.length):
            self.weights.append(float(rate))
        return self.weights
