#!/usr/bin/env python

from collections import Counter
import numpy as np
import simulate
import evaluate

#-------------------------------------------------------------------------------


#set rates for transcript
rates = {}
rates["TAT"] = 1.0
rates["TAC"] = 0.5
sequences = list(rates.keys())
rates= list(rates.values())

#generate and save randomized transcripts
gene_rates = []
for i in range(0,120):
    r = np.random.randint(0,2)
    gene_rates.append(float(rates[r]))
gene_rates = np.repeat(gene_rates,3)
gene_length = len(gene_rates)
rate_counts = Counter(gene_rates)

with open('transcripts.tsv','w') as summary:
    summary.writelines(str(generation)+ '\t' + str(gene_rates)+'\n')
summary.close()
print(gene_rates)

#run pinetree
if __name__ == "__main__":
    simulate.execute("script_function_test",gene_rates,gene_length)

#evaluate protein production
simulation = open("fast_transcript_3_counts.tsv", 'r')
proteins = []
time = []
feature = []
lines = simulation.readlines()
linecounter = 0
for line in lines:
    if linecounter == 0:
        print("header found")
    else:
        data = line.strip().split('\t')
        proteins.append(float(data[2]))
        time.append(float(data[0]))
        feature.append(str(data[1]))
    linecounter += 1
simulation.close()

#filter table for protein feature
ribosome_counts = []
protein_counts = []
times = []
for i in range(len(feature)):
    if feature[i] =='__ribosome':
        ribosome_counts.append(proteins[i])
        times.append(time[i])
    if feature[i] == 'proteinP':
        protein_counts.append(proteins[i])

#calculate steady state using evaluate.py
list_len = len(times)
slopes = []
evaluate.slope(times,slopes,ribosome_counts)
evaluate.steady_state(times,slopes,ribosome_counts,protein_counts)

def safe_calc(exponent):
  if exponent > 700:
    print("system maxed")
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))

# accelerated simulation functions
def calc_prob_scores(stab_mut, stab_org,N):

xi = stab_org
xj = stab_mut

  if xj >= xi:
    return(float(1.0))
  else:
exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))
