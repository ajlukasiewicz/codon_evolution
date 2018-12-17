#!/usr/bin/env python

import pinetree as pt
import numpy as np
from collections import Counter

rates = {}
rates["TAT"] = 1.0
rates["TAC"] = 0.5
sequences = list(rates.keys())
rates= list(rates.values())
#print(sequences)

gene_rates = []
for i in range(0,120):
    r = np.random.randint(0,2)
    gene_rates.append(float(rates[r]))
# for i in range(0,120):
#     gene_rates.append(1.0)
gene_rates = sorted(gene_rates*3)
gene_length = len(gene_rates)

print(gene_rates)
rate_counts = Counter(gene_rates)
print(rate_counts)

def execute(output):
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_ribosome(copy_number=100,speed=30,footprint=10)
    transcript = pt.Transcript(name="transcript",length=gene_length)
    transcript.add_gene(name="proteinP", start=30, stop=360,rbs_start=(30 - 15), rbs_stop=30, rbs_strength=1e7)
    transcript.add_weights(weights=gene_rates)
    sim.register_transcript(transcript)
    sim.simulate(time_limit=200, time_step=1, output=output + "_counts.tsv")

if __name__ == "__main__":
    execute("fast_transcript_3")
