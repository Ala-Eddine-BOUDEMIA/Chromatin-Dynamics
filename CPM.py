import pandas as pd

import Config 

def generate_cpm(raw, normalized):

    raw_counts = pd.read_csv(raw, 
        header = 0, index_col = 0, sep = '\t')
	
    # Sum of counts per sample
    counts_per_sample = raw_counts.sum(axis = 0)
    print(counts_per_sample)

    # cpm
    total = counts_per_sample.div(1e6)
    cpm = raw_counts.loc[:,:].div(total) 

    print(cpm.sum(axis = 0))
    cpm.to_csv(str(normalized), sep = "\t")

if __name__ == '__main__':
    generate_cpm(
        raw = Config.args.bf,
        normalized = Config.args.norm)