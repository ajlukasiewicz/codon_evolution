#!/usr/bin/env python

# Author: Alexandra Lukasiewicz
# Date: 1/10/2019

import argparse
from decimal import Decimal, ROUND_HALF_EVEN
import numpy as np
from time import time

from transcripts import Transcript
import ptsim
import fit_eval


# -----------------------------------------------------------------------------


# Introduce mutation in transcript
def mutate(transcript):
    rates = [1.0, 100.0]
    n = np.random.randint(0, len(transcript))
    m = np.random.randint(0, 2)
    transcript[n] = rates[m]
    return transcript


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
def write_transcript_data(transcript_data, transcript, slow_rt, fast_rt, production_rate):
    transcript_data.append(str(gen) + '_' + str(population) + ','
                           + str(slowcodon_percent) + ','
                           + str(fastcodon_percent) + ','
                           + str(rate) + ','
                           + ','.join(map(str,transcript)) + '\n')
    return transcript_data

def main():
    parser = argparse.ArgumentParser(description='run simulation for x generations')
    parser.add_argument(
         '-g',
         action='store',
         dest='g',
         required=True,
         type=int,
         help="number of generations to evolve for",
         )
    options = parser.parse_args()

    max_generations = options.g
    
    t0=time()

    i = 0
    gen = 0
    N = 2   # population size

    # Write data headers
    transcript_data = ['slow_rt,fast_rt,prod_rt,transcript \n']

    rates = [1,100]

    for n in range(1,101):
        rates[0] = n
        g1 = Transcript(120,1,rates)
        g1_weights = g1.random_codons()

        # Create and evaluate original transcript
        ptsim.simulate('_'.join(map(str,rates)), 'g1', g1_weights)
        g1_fitness = Decimal(fit_eval.main("generation_" + '_'.join(map(str,rates)) + '_' + 'g1' + '_counts.tsv', 'proteinP')).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)
        transcript_data.append(str(rates[0]) + ',' + str(rates[1]) + ','
                           + str(g1_fitness) + ','
                           + ','.join(map(str,g1_weights)) + '\n')

    with open('transcript_stats.csv', 'w') as f:
        f.writelines(transcript_data)
        f.close()
    
    t1 = time()
    print('program takes %f' %(t1-t0))

if __name__ == '__main__':
    main()