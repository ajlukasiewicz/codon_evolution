#!/usr/bin/env python

# Author: Alexandra Lukasiewicz
# Date: 1/10/2019

import argparse
from decimal import Decimal, ROUND_HALF_EVEN
import numpy as np
from time import time
import sys

from transcripts import Transcript
import pooled_transcripts as pt
import fit_eval

# -----------------------------------------------------------------------------


# Introduce mutation in transcript
def mutate(transcript,min,max):
    rates = [min, max]
    n = np.random.randint(0, len(transcript))
    m = np.random.randint(0, 2)
    transcript[n] = rates[m]
    return transcript, n, rates[m]


# Accelerated simulation functions
def safe_calc(exponent):
    if exponent > 700:
        print("system maxed")
        return(sys.float_info.max)
    else:
        return(Decimal(exponent).exp())


def calc_prob_scores(stab_mut, stab_org, N):
    xi = stab_org
    xj = stab_mut

    if xj >= xi:
        return(1.0)
    else:
        exponent = -2 * Decimal(N) * (xi - xj)
        return(safe_calc(exponent))


# Write transcript metrics to tsv output
def write_transcript_data(transcript_data, transcript, gen, population, rate, mutation_loc, mutation, fast_speed, slow_speed):
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

def write_speed_locations(gen,transcript,location_summary, min):
    rbs = transcript[0:11]
    orf = transcript[11:71]
    term = transcript[71:120]
    location_summary.append(str(gen) + ',' 
                            + str(rbs.count(min)) + ',' 
                            + str(orf.count(min)) + ',' 
                            + str(term.count(min)) + '\n')
    return location_summary

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
         '-r',
         action='store',
         dest='r',
         required=False,
         default=0.5,
         type=float,
         help="Set the slow rate for the system, default = 0.5",
         )
    parser.add_argument(
         '-m',
         action='store',
         dest='m',
         required=False,
         default=1.0,
         type=float,
         help="Set the fast rate for the system, default = 1.0",
         )
      
    options = parser.parse_args()
    max_generations = options.g
    t0 = time()

    i = 0
    gen = 0
    N = 10   # population size

    # Write data headers
    transcript_data = ['Generation,Population,Slow_Count,Fast_Count,Mut_Loc,Prod_rt,Min \n']
    location_summary = ['Generation,RBS,ORF,Terminator \n']


    rates = [options.r, options.m]

    while i <= max_generations:

        # Initialize model with random transcript
        if gen == 0:
            popA = Transcript(120, gen, rates)
            popA_weights = popA.random_codons()

        # Mutate original transcript and evaluate performance
        popB_weights, mutation_loc, mutation = mutate(popA_weights[:],options.r,options.m)

        # Create and evaluate original transcript
        pt.simulate(gen, popA_weights, popB_weights, rates)
        popA_fitness = Decimal(fit_eval.main("../data/generation_" + str(gen) + '_' + str(rates[0]) + '_' + str(rates[1]) + '_counts.tsv', 'proteinA')).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)
        popB_fitness = Decimal(fit_eval.main("../data/generation_" + str(gen)+ '_' +  str(rates[0]) + '_' + str(rates[1]) + '_counts.tsv', 'proteinB')).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)
        
        p = Decimal(np.random.uniform()).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)
        print(p)
        # Compare production rates and calculate probability of mutation acceptance
        probability = Decimal(calc_prob_scores(popB_fitness, popA_fitness, N)).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)

        if probability == Decimal('1.000') or p <= probability:
            write_transcript_data(transcript_data, popB_weights, gen, 'popB', popB_fitness, mutation_loc, mutation, options.m, options.r)
            popA_weights = popB_weights[:]
            write_speed_locations(gen,popB_weights,location_summary, options.r)

        else:
            write_transcript_data(transcript_data, popA_weights, gen, 'popA', popA_fitness, 0, 0, options.m, options.r)
            popA_weights = popA_weights
            write_speed_locations(gen,popA_weights,location_summary, options.r)

        if gen == max_generations:
            final_fit = str(rates[0]) + ',' + str(rates[1]) + ',' + str(popA_fitness) + '\n'

        gen += 1
        i += 1

    with open('../data/transcript_stats_' + str(rates[0]) + '_' +  str(rates[1]) + '.csv', 'w') as f:
        f.writelines(transcript_data)
        f.close()
    
    with open('../data/final_fitness_' + str(rates[0]) + '_' +  str(rates[1]) + '.csv', 'w') as f:
        f.writelines(final_fit)
        f.close()
    
    with open('../data/location_stats_'+ str(rates[0]) + '_' +  str(rates[1]) + '.csv','w') as f:
        f.writelines(location_summary)
        f.close()

    t1 = time()
    print('program takes %f' % (t1-t0))


if __name__ == '__main__':
    main()
