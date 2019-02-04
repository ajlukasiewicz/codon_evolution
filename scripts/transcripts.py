#!/usr/bin/env python
'''
Define codon translation rates for use in pinetree simulation
Author: Alexandra Lukasiewicz
'''
import numpy as np


# ------------------------------------------------------------------------------
class Transcript:

    def __init__(self, length, generation, rates):
        self.length = length
        self.weights = []
        self.generation = generation
        self.rates = rates
 
    def random_codons(self):
        for i in range(0, self.length):
            r = np.random.randint(0, 2)
            self.weights.append(self.rates[r])
        return self.weights

    def defined_codons(self, rate):
        for i in range(0, self.length):
            self.weights.append(float(rate))
        return self.weights
    
#    def vary_rate(self, max, min):
#        for n in range(min, max):
#            self.rates[0] = n
#            self.weights = self.weights.random_codons()
#        return self.weights