#!/usr/bin/env python

# Author: Alexandra Lukasiewicz
# Date: 1/10/2019

import argparse
from decimal import Decimal, ROUND_HALF_EVEN
import numpy as np
import math

from transcripts import Transcript
import ptsim
import fit_eval

# -----------------------------------------------------------------------------


# Introduce mutation in transcript
def mutate(transcript):
    rates = [0.5, 1.0]
    n = np.random.randint(0, len(transcript))
    m = np.random.randint(0, 2)
    transcript.insert(n, rates[m])
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
def write_transcript_data(transcript_data, transcript, gen, population, rate):
    slowcodons = transcript.count(0.5)
    fastcodons = transcript.count(1.0)
    slowcodon_percent = (slowcodons / len(transcript))*100
    fastcodon_percent = (fastcodons / len(transcript))*100
    transcript_data.append(str(gen) + '_' + str(population) + ','
                           + str(slowcodon_percent) + ','
                           + str(fastcodon_percent) + ','
                           + str(rate)
                           + '\n')
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

    i = 0
    gen = 0
    N = 2   # population size

    # Write data headers
    transcript_data = ['Generation,Pct_Slow,Pct_Fast,Rate \n']

    while i <= max_generations:

        # Initialize model with random transcript
        if gen == 0:
            popA = Transcript(120, gen)
            popA_weights = popA.random_codons()

        # Create and evaluate original transcript
        ptsim.simulate(gen, 'popA', popA_weights)
        popA_fitness = Decimal(fit_eval.main("generation_" + str(gen)+ '_' + 'popA' + '_counts.tsv', 'proteinP')).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)
        write_transcript_data(transcript_data, popA_weights, gen, 'popA', popA_fitness)
        print(popA_fitness)

        # Mutate original transcript and evaluate performance
        popB_weights = mutate(popA_weights[:])
        ptsim.simulate(gen, 'popB',popB_weights)
        popB_fitness = Decimal(fit_eval.main("generation_" + str(gen)+ '_' + 'popB' + '_counts.tsv', 'proteinP')).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)
        write_transcript_data(transcript_data, popB_weights, gen, 'popB', popB_fitness)
        print(popB_fitness)

        # Compare production rates and calculate probability of mutation acceptance
        probability = Decimal(calc_prob_scores(popB_fitness, popA_fitness, N)).quantize(Decimal('.001'), rounding = ROUND_HALF_EVEN)
        print(probability)

        if probability >= Decimal('1.000'):
            popA_weights = popB_weights[:]

        else:
            popA_weights = popA_weights

        gen += 1
        i += 1

    with open('transcript_stats.csv', 'w') as f:
        f.writelines(transcript_data)
        f.close()


if __name__ == '__main__':
    main()
