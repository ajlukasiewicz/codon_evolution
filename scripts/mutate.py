#!/usr/bin/env python

import numpy as np
import transcripts


# mutate
def mutate(transcript):
    n = np.random.randint(0,len(gene))
    m = randint(0,2)
    transcript.insert(n,rates[m])


# accelerated simulation functions

def safe_calc(exponent):
    if exponent > 700:
        print("system maxed")
        return(sys.float_info.max)
    else:
        return(math.exp(exponent))

def calc_prob_scores(stab_mut, stab_org,N):

xi = stab_org
xj = stab_mut

  if xj >= xi:
    return(float(1.0))
  else:
exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))
