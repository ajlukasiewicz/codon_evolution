#!/usr/bin/env python

import numpy as np

import transcripts
from transcripts import Transcript
import ptsim
import fit_eval

#-------------------------------------------------------------------------------

# Mutate codon speed at a single random location in transcript
def mutate(transcript):
    rates = [0.5,1.0]
    n = np.random.randint(0,len(transcript))
    m = np.random.randint(0,2)
    transcript.insert(n,rates[m])

# accelerated simulation functions
# def safe_calc(exponent):
#     if exponent > 700:
#         print("system maxed")
#         return(sys.float_info.max)
#     else:
#         return(math.exp(exponent))
#
# def calc_prob_scores(stab_mut, stab_org,N):
#
# xi = stab_org
# xj = stab_mut
#
#   if xj >= xi:
#     return(float(1.0))
#   else:
# exponent = -2 * float(N) * (xi - xj)
#     return(safe_calc(exponent))

# Write transcript metrics to tsv output
def write_transcript_data(transcript,gen,rate):
    transcript_data = ['Generation','% Slow','% Fast','rate','transcript \n']
    slowcodons = ((1 for codon in transcript if codon == 0.5)
                   / len(transcript))*100
    fastcodons = ((1 for codon in trasnscript if codon == 1.0)
                   / len(transcript))*100
    transcript_data.append(str(gen) + ','
                           + str(slowcodons) + ','
                           + str(fastcodons) + ','
                           + str(rate) + ','
                           + str(transcript) + '\n')

    with open('transcript_stats.csv', 'w') as f:
        f.writelines(transcript_data)
    f.close()

def main():
    max_generations = 2
    i = 0
    gen=0
    while i <= max_generations:
        gen+=1
        genA = Transcript(120,gen)
        genA_weights = genA.random_codons()
        genB = mutate(genA_weights)
        ptsim.simulate(gen,genA_weights)
        fit_eval.main("generation_" + str(gen)+ '_counts.tsv', 'proteinP')

        i+=1
if __name__ == '__main__':
    main()
