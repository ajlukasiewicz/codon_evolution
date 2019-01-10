#!/usr/bin/env python3


import evaluate as ev

#-------------------------------------------------------------------------------

def main(simulation_output,protein_name):
    
    g1 = ev.EvalSteadyState()
    g1.importTable(simulation_output,protein_name)
    g1.ribo_slope(g1.time,g1.ribosome_counts)
    g1.steady_state(g1.time,g1.slopes,g1.ribosome_counts,g1.protein_counts)
    g1.linreg(g1.x,g1.y)
    fitness = g1.fitness
    return fitness

if __name__ == "__main__":
     main()
