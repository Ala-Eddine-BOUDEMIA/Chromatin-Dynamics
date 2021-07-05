import pandas as pd 
from numpy.random import default_rng

from pathlib import Path

import Config

def generate_data(
	cv, meta, rand, top88, cv_list, nonRcv,
	top1000, tissue_counts, normalized_counts,
	nonRcv_list):
	
	counts = pd.read_csv(normalized_counts,
		header = 0, index_col = 0, sep = "\t")

	chaperones_variants = pd.read_csv(cv_list,
		header = 0, index_col = 0, sep = ";")

	chaperone_nonRv = pd.read_csv(nonRcv_list,
		header = 0, index_col = 0, sep = "\t")

	metadata = pd.read_csv(meta,
		header = 0, index_col = 0, sep = "\t")
	
	# Generate the top 88 expressed genes
	# Generate the top1000 expressed genes
	
	counts["total"] = counts.sum(axis = 0)
	counts = counts.sort_values("total")
	counts.pop("total")
	
	top88_g = counts.iloc[:88, :]
	top88_g.to_csv(top88, sep = "\t")
	
	top1000_g = counts.iloc[:1000, :]
	top1000_g.to_csv(top1000, sep = "\t")
	
	# Generate chaperones and variants dataframe
	# Generate chaperones and non replicative variants dataframe
	for i in counts.index.to_list():
		counts.rename(index = {i: i.split('.')[0]}, inplace = True)
	
	cv_df = pd.DataFrame(columns = counts.columns)
	nonRv_df = pd.DataFrame(columns = counts.columns)
	
	for i in counts.index.to_list():
		
		for j in chaperones_variants.index.to_list():
			if i.strip() == j.strip():
				cv_df = cv_df.append(counts.loc[i])
	
		for k in chaperone_nonRv.index.to_list():
			if i.strip() == k.strip():
				nonRv_df = nonRv_df.append(counts.loc[i])

	cv_df.to_csv(str(cv), sep = "\t")
	nonRv_df.to_csv(str(nonRcv), sep = "\t")
	
	# Generate random sets
	for c in range(10):
		rng = default_rng()
		r = rng.choice(len(counts), size = 88, replace = False)
		df_random = pd.DataFrame(columns = counts.columns)
		
		for i in r:
			df_random = df_random.append(counts.iloc[i])

		df_random.to_csv(rand.joinpath(str(c) + ".tsv"), sep = "\t")

	# Generate counts by tissue
	counts = counts.T
	counts = counts.join(metadata["smts"])
	tissue_types = pd.unique(metadata["smts"])

	for t in tissue_types:
		df = pd.DataFrame(columns = counts.columns)
		df = df.append(counts[counts["smts"] == t])
		df.pop("smts")
		df.to_csv(tissue_counts.joinpath(t + ".tsv"), sep = '\t')

if __name__ == '__main__':

	generate_data(
		cv = Config.args.cv,
		meta = Config.args.meta,
		rand = Config.args.rand,
		top88 = Config.args.top88,
		cv_list = Config.args.list,
		nonRcv = Config.args.nonRcv,
		top1000 = Config.args.top1000,
		tissue_counts = Config.args.tissue,
		normalized_counts = Config.args.norm,
		nonRcv_list = Config.args.nonReplicative)