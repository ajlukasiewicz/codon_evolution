#!/usr/bin/env python3

import argparse

import evaluate as ev

#-------------------------------------------------------------------------------

def main(simulation_output,protein_name):
    # parser = argparse.ArgumentParser(description='simulation counts table')
    # parser.add_argument(
    #      '-i',
    #      action='store',
    #      dest='i',
    #      required=True,
    #      type=str,
    #      help="input .tsv file from pinetree",
    #      )
    # parser.add_argument(
    #      '-f',
    #      action='store',
    #      dest='f',
    #      required=True,
    #      type=str,
    #      help="name of protein in system",
    #      )
    # options = parser.parse_args()

    g1 = ev.EvalSteadyState()
    g1.importTable(simulation_output,protein_name)
    g1.ribo_slope(g1.time,g1.ribosome_counts)
    g1.steady_state(g1.time,g1.slopes,g1.ribosome_counts,g1.protein_counts)
    g1.linreg(g1.x,g1.y)

if __name__ == "__main__":
     main()
