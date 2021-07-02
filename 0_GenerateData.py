import pandas as pd 
from numpy.random import default_rng

from pathlib import Path

import Config

def generate_data(
	normalized_counts, top1000, top81, 
	rand, cv, nonRcv, cv_list, nonRcv_list):
	
	counts = pd.read_csv(normalized_counts,
		header = 0, index_col = 0, sep = "\t")

	chaperones_variants = pd.read_csv(cv_list,
		header = 0, index_col = 0, sep = ";")

	chaperone_nonRv = pd.read_csv(nonRcv_list,
		header = 0, index_col = 0, sep = "\t")
	
	# Generate the top 81 expressed genes
	# Generate the top1000 expressed genes
	
	counts["total"] = counts.sum(axis = 0)
	counts = counts.sort_values("total")
	counts.pop("total")
	
	top81_g = counts.iloc[:81, :]
	top81_g.to_csv(top81, sep = "\t")
	
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
		r = rng.choice(len(counts), size = 81, replace = False)
		df_random = pd.DataFrame(columns = counts.columns)
		
		for i in r:
			df_random = df_random.append(counts.iloc[i])

		df_random.to_csv(rand.joinpath(str(c) + ".tsv"), sep = "\t")

if __name__ == '__main__':
	generate_data(
		normalized_counts = Config.args.norm,
		top1000 = Config.args.top1000,
		top81 = Config.args.top81,
		rand = Config.args.rand,
		cv = Config.args.cv,
		nonRcv = Config.args.nonRcv,
		cv_list = Config.args.list,
		nonRcv_list = Config.args.nonReplicative)