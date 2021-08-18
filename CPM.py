import pandas as pd

import Config 

def generate_cpm(meta, raw, normalized):

    raw_counts = pd.read_csv(raw, 
        header = 0, index_col = 0, sep = '\t')
	
    # Sum of counts per sample
    counts_per_sample = raw_counts.sum(axis = 0)
    print(counts_per_sample)

    # cpm
    total = counts_per_sample.div(1e6)
    cpm = raw_counts.loc[:,:].div(total) 

    print(cpm.sum(axis = 0))
    cpm['total'] = cpm.sum(axis = 1)
    cpm = cpm.drop(cpm[cpm["total"] <= 1].index)
    #cpm = cpm.drop(cpm[cpm["total"] > 4000].index)
    cpm.pop('total')

    if Config.args.dataset == 'TCGA':
        metadata = pd.read_csv(meta, 
            header = 0, index_col = 0, sep = '\t')
        cpm = cpm.T
        cpm = cpm.join(metadata['gdc_cases.samples.sample_type'])
        cpm = cpm[cpm['gdc_cases.samples.sample_type'] == 'Solid Tissue Normal']
        cpm.pop('gdc_cases.samples.sample_type')
        cpm = cpm.T

    cpm.to_csv(str(normalized), sep = "\t", 
        float_format='%.3f')

if __name__ == '__main__':
    generate_cpm(
        meta = Config.args.meta,
        raw = Config.args.bf,
        normalized = Config.args.norm)