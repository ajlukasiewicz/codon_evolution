import pandas as pd
import numpy as np

def mega_import(foldername):
    ribo_production = pd.read_csv("../data/" + str(foldername) + "/transcript_stats_0.5_1.0.csv")
    return ribo_production

folder_names = []
ribo_production_all = {}
n = 0 
for i in range(5,20):
    folder_names.append("prod_rates_" + str(i))
    ribo_production_all[folder_names[n]] = mega_import("ribosome_num_" + str(i))
    n += 1

print(ribo_production_all)
#for k,v in ribo_production_all.items():
