#!/usr/bin/env python
'''
Simulate protein production using pinetree
Author: Alexandra Lukasiewicz
'''
import numpy as np

import pinetree as pt

# import transcripts
# from transcripts import Transcript

# ------------------------------------------------------------------------------


def execute(outfile, gene_rates, gene_length):
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_ribosome(copy_number=100, speed=30, footprint=10)
    transcript = pt.Transcript(name="transcript", length=gene_length)
    transcript.add_gene(name="proteinP", start=30, stop=360, rbs_start=(30 - 15), rbs_stop=30, rbs_strength=1e7)
    transcript.add_weights(weights=gene_rates)
    sim.register_transcript(transcript)
    sim.simulate(time_limit=100, time_step=1, output=outfile + "_counts.tsv")


def simulate(gen, pop, transcript_weights):
    transcript_weights = np.repeat(transcript_weights, 3)
    transcript_length = len(transcript_weights)
    execute("generation_" + str(gen) + '_' + str(pop), transcript_weights, transcript_length)

if __name__ == "__main__":
    main()
