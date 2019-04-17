#!/usr/bin/env python
'''
Simulate protein production using pinetree
Author: Alexandra Lukasiewicz
'''
import numpy as np
from transcripts import *
import pinetree as pt

# import transcripts
# from transcripts import Transcript

# ------------------------------------------------------------------------------


def execute(outfile, gene_rates_A, gene_rates_B, gene_length):
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_ribosome(copy_number=5, speed=30, footprint=10)
    transcript_A = pt.Transcript(name="transcript_A", length=gene_length)
    transcript_A.add_gene(name="proteinA", start=30, stop=360, rbs_start=(30 - 15), rbs_stop=30, rbs_strength=1e7)
    transcript_A.add_weights(weights=gene_rates_A)
    transcript_B = pt.Transcript(name="transcript_B", length=gene_length)
    transcript_B.add_gene(name="proteinB", start=30, stop=360, rbs_start=(30 - 15), rbs_stop=30, rbs_strength=1e7)
    transcript_B.add_weights(weights=gene_rates_B)
    sim.register_transcript(transcript_A)
    sim.register_transcript(transcript_B)
    sim.simulate(time_limit=100, time_step=1, output=outfile + "_counts.tsv")


def simulate(gen,transcript_weights_A, transcript_weights_B,rates):
    transcript_A = np.repeat(transcript_weights_A, 3)
    transcript_B = np.repeat(transcript_weights_B, 3)
    transcript_length = len(transcript_A)
    execute("../data/generation_" + str(gen) + '_' + str(rates[0]) + '_' + str(rates[1]), transcript_A, transcript_B, transcript_length)


#def main():
#    rates = [0.5, 1.0]
#    poolA = Transcript(120,0,rates)
#    poolA_weights = poolA.random_codons()
#    poolB_weights, mutation_loc, mutation = mutate(poolA_weights[:],0.5,1.0)
#    #simulate the pooled system here
#    simulate(poolA_weights,poolB_weights) 

if __name__ == "__main__":
    main()