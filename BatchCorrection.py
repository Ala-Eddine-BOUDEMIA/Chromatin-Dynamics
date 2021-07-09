import pandas as pd 
import numpy as np

from combat.pycombat import pycombat

import Config

def batch_correction(meta, filtered_counts):
    
    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = "\t")

    counts = pd.read_csv(filtered_counts, header = 0, index_col = 0, sep = "\t")
    counts = pd.DataFrame(np.log2(counts + 1))

    batch1 = metadata["smnabtch"].astype(str)
    batch2 = metadata["smgebtch"].astype(str)

    counts_corrected1 = pycombat(counts, batch1)
    counts_corrected1 = counts_corrected1.fillna(0)

    counts_corrected = pycombat(counts_corrected1, batch2)
    counts_corrected = counts_corrected.fillna(0)

    counts_corrected.to_csv("counts_filtered_corrected.tsv", sep = "\t")

if __name__ == '__main__':	

    batch_correction(
        meta = Config.args.meta,
        filtered_counts = Config.args.af)