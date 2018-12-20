
import numpy as np
import pinetree as pt
import transcripts
from transcripts import Transcript


def execute(outfile,gene_rates,gene_length):
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_ribosome(copy_number=100,speed=30,footprint=10)
    transcript = pt.Transcript(name="transcript",length=gene_length)
    transcript.add_gene(name="proteinP", start=30, stop=360,rbs_start=(30 - 15), rbs_stop=30, rbs_strength=1e7)
    transcript.add_weights(weights=gene_rates)
    print("added weights")
    sim.register_transcript(transcript)
    print("registered transcript")
    sim.simulate(time_limit=100, time_step=1, output=outfile + "_counts.tsv")


def main():
    test_transcript = Transcript(120)
    print(test_transcript)
    test_transcript_weights = test_transcript.random_codons()
    test_transcript_weights = np.repeat(test_transcript_weights,3)
    print(test_transcript_weights)
    transcript_length = len(test_transcript_weights )#include this in simulation
    execute("test_counts" , test_transcript_weights, transcript_length)

if __name__ == "__main__":
    main()
