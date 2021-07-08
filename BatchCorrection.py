import pandas as pd

from combat.pycombat import pycombat

import Config

metadata = pd.read_csv(Config.args.meta, 
	header = 0, index_col = 0, sep = "\t")

counts = pd.read_csv(Config.args.af,
	header = 0, index_col = 0, sep = "\t")

counts = counts + 0.1 

batch = metadata["smnabtch"]

df_corrected = pycombat(counts, batch)

df_corrected.to_csv("Corrected_.tsv", sep = "\t")