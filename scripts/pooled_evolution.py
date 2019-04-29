#!/usr/bin/env python

# Author: Alexandra Lukasiewicz
# Date: 1/10/2019

import argparse
import numpy as np
import math
import subprocess
from statistics import mean
from time import time
import sys

from transcripts import Transcript
import pooled_transcripts as pt
import fit_eval

# -----------------------------------------------------------------------------


# Introduce mutation in transcript
def mutate(transcript, min, max):
    transcript2 = transcript.copy()
    rates = [min, max]
    n = np.random.randint(0, len(transcript2))
    m = np.random.randint(0, 2)
    transcript2[n] = rates[m]
    return transcript2, n, rates[m]


# Accelerated simulation functions
def safe_calc(exponent):
    if exponent > 700:
        print("system maxed")
        return(sys.float_info.max)
    else:
        return(math.exp(exponent))


def calc_prob_scores(stab_mut, stab_org, N):
    xi = stab_org
    xj = stab_mut

    if xj >= xi:
        return(1.0)
    else:
        exponent = -2 * N * (xi - xj)
        return(safe_calc(exponent))

def replicate_simulations(gen, popA_weights, popB_weights, rates, ribosomes, speed, outfile):
    popA_fit_reps = []
    popB_fit_reps = []
    for n in range(0,4):
        pt.simulate(gen, popA_weights, popB_weights, rates, ribosomes, speed, outfile, n)
        popA_fit_reps.append(fit_eval.main("../data/" + str(outfile) + "/generation_" + str(gen) + '_' + 'replicate_' + str(n) + '_' + str(rates[0]) + '_' + str(rates[1]) + '_counts.tsv', 'proteinA'))
        popB_fit_reps.append(fit_eval.main("../data/" + str(outfile) + "/generation_"  + str(gen)+ '_'  + 'replicate_' + str(n) + '_' + str(rates[0]) + '_' + str(rates[1]) + '_counts.tsv', 'proteinB'))
    popA_finess = mean(popA_fit_reps)
    popB_fitness = mean(popB_fit_reps)
    return popA_finess, popB_fitness

# Write transcript metrics to tsv output
def write_transcript_data(transcript_data, transcript, gen, population, rate, mutation_loc, mutation, slow_speed, fast_speed):
    fastcodons = transcript.count(fast_speed)
    slowcodons = len(transcript) - fastcodons
    transcript_data.append(str(gen) + ',' + str(population) + ','
                           + str(slowcodons) + ','
                           + str(fastcodons) + ','
                           + str(mutation_loc) + ','
                           + str(mutation) + ','
                           + str(rate) +  ','
                           + str(slow_speed) + '\n')
    return transcript_data


def main():
    parser = argparse.ArgumentParser(description='run simulation for x generations')
    parser.add_argument(
         '-g',
         action='store',
         dest='g',
         required=True,
         type=int,
         help="Set the number of generations to evolve for",
         )
    parser.add_argument(
         '-s',
         action='store',
         dest='s',
         required=False,
         default=0.5,
         type=float,
         help="Set the slow rate for the system, default = 0.5",
         )
    parser.add_argument(
         '-f',
         action='store',
         dest='f',
         required=False,
         default=1.0,
         type=float,
         help="Set the fast rate for the system, default = 1.0",
         )
    parser.add_argument(
         '-r',
         action='store',
         dest='r',
         required=False,
         default=5,
         type=int,
         help="Set number of ribosomes in simulation, default = 5",
         )
    parser.add_argument(
         '-sp',
         action='store',
         dest='sp',
         required=False,
         default=30,
         type=int,
         help="Set the ribosome speed in simulation, default = 30",
         )
    parser.add_argument(
         '-o',
         action='store',
         dest='o',
         required=True,
         default=30,
         type=str,
         help="output directory name",
         )
      
    options = parser.parse_args()
    max_generations = options.g
    t0 = time()

    subprocess.call("mkdir ../data/" + str(options.o), shell = True)

    i = 0
    gen = 0
    N = 10000   # population size

    # Write data headers
    transcript_data = ['Generation,Population,Slow_Count,Fast_Count,Mut_Loc,Mutation,Prod_rt,Min \n']
    outfasta = []

    rates = [options.s, options.f]
    
    # Initialize model with random transcript
    if gen == 0:
            popA = Transcript(120, gen, rates)
            dist = [0.5,0.5]
            popA_weights = popA.random_with_dist(rates,dist)
        # Mutate original transcript and evaluate performance
            popB_weights, mutation_loc, mutation = mutate(popA_weights, options.s, options.f)
            print(popA_weights,popB_weights)


    while i <= max_generations:
        # Create and evaluate original transcript
        popA_fitness, popB_fitness = replicate_simulations(gen, popA_weights, popB_weights, rates, options.r, options.sp , options.o)
        print(popA_fitness, popB_fitness)
        
        # Compare production rates and calculate probability of mutation acceptance
        p = np.random.uniform()
        probability = calc_prob_scores(popB_fitness, popA_fitness, N)
        print(p,probability)

        if p <= probability:
            write_transcript_data(transcript_data, popB_weights, gen, 'popB', popB_fitness, mutation_loc, mutation, options.s, options.f)
            outfasta.append(">Generation " + str(gen) + "\n" + str(popB_weights) + "\n" )
            popA_weights = popB_weights
            popB_weights, mutation_loc, mutation = mutate(popA_weights, options.s, options.f)

            
        else:
            write_transcript_data(transcript_data, popA_weights, gen, 'popA', popA_fitness, 0, 0, options.s, options.f)
            outfasta.append(">Generation " + str(gen) + "\n" + str(popA_weights) + "\n" )
            popA_weights = popA_weights
            popB_weights, mutation_loc, mutation = mutate(popA_weights, options.s, options.f)

            
        print(popA_weights,popB_weights)
        if gen == max_generations:
            final_fit = str(rates[0]) + ',' + str(rates[1]) + ',' + str(popA_fitness) + '\n'

        gen += 1
        i += 1

    with open('../data/' + str(options.o) + '/transcript_stats_' + str(rates[0]) + '_' +  str(rates[1]) + '.csv', 'w') as f:
        f.writelines(transcript_data)
        f.close()
    
    with open('../data/' + str(options.o) + '/transcripts_' + str(rates[0]) + '_' +  str(rates[1]) + '.fasta', 'w') as f:
        f.writelines(outfasta)
        f.close()

    t1 = time()
    print('program takes %f' % (t1-t0))


if __name__ == '__main__':
    main()
