
import argparse
from decimal import Decimal, ROUND_HALF_EVEN
import numpy as np
from time import time
import sys

from transcripts import Transcript
import ptsim
import evaluate as ev

gen = 1
ratesA = [1,1]
ratesB = [5,5] 

popA = Transcript(120, gen, ratesA)
popA_weights = popA.defined_codons(1)
popB = Transcript(120, gen, ratesB)
popB_weights = popB.defined_codons(5)

ptsim.simulate(gen,'popA',popA_weights, ratesA)
ptsim.simulate(gen,'popB',popB_weights, ratesB)

print(popA_weights)
print(popB_weights)

#evaluate fitness
popA = ev.EvalSteadyState()
popA.importTable('../data/generation_1_popA_1_1_counts.tsv', 'proteinP')
popA.ribo_slope(popA.time, popA.ribosome_counts)
popA.steady_state(popA.time, popA.slopes, popA.ribosome_counts, popA.protein_counts)
popA.linreg(popA.x, popA.y)
popA_fitness = popA.fitness

print(popA_fitness)

popB = ev.EvalSteadyState()
popB.importTable('../data/generation_1_popB_5_5_counts.tsv', 'proteinP')
popB.ribo_slope(popB.time, popB.ribosome_counts)
popB.steady_state(popB.time, popB.slopes, popB.ribosome_counts, popB.protein_counts)
popB.linreg(popB.x, popB.y)
popB_fitness = popB.fitness

print(popB_fitness)