import pandas as pd 

from combat.pycombat import pycombat

import Config

def batch_correction(meta, filtered_counts):
    
    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = "\t")

    counts = pd.read_csv(filtered_counts, header = 0, index_col = 0, sep = "\t")
    counts = counts + 1
    print(counts.head())

    batch1 = metadata["smnabtch"]
    batch2 = metadata["smgebtch"]

    counts_corrected_batch1 = pycombat(counts, batch1)
    print(counts_corrected_batch1.head())
    counts_corrected = pycombat(counts_corrected_batch1, batch2)
    print(counts_corrected.head())

    counts_corrected.to_csv("counts_filtered_corrected.tsv", sep = "\t")

if __name__ == '__main__':	

    batch_correction(
        meta = Config.args.meta,
        filtered_counts = Config.args.af)