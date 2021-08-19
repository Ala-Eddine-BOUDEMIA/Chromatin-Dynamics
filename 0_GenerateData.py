import pandas as pd 
from numpy.random import default_rng

from pathlib import Path

import Tools
import Config

def generate_data(
	full, meta, rand, top100, cv_list, 
	nonRcv, top1000, normal, tissue_counts, 
	normalized_counts, nonRcv_list, wo_bbbpst_tissues):
	
	counts = pd.read_csv(normalized_counts,
		header = 0, index_col = 0, sep = "\t")

	chaperones_variants = pd.read_csv(cv_list,
		header = 0, index_col = 0, sep = ";")

	chaperone_nonRv = pd.read_csv(nonRcv_list,
		header = 0, index_col = 0, sep = "\t")

	metadata = pd.read_csv(meta,
		header = 0, index_col = 0, sep = "\t")
	
	# Generate the top 125 expressed genes
	# Generate the top1000 expressed genes
	Tools.create_folder(("/").join(str(top100).split("/")[:-1]))
	Tools.create_folder(("/").join(str(top1000).split("/")[:-1]))

	counts["total"] = counts.sum(axis = 1)
	counts = counts.sort_values("total", ascending = False)
	counts.pop("total")
	
	top100_g = counts.iloc[:125, :]
	top100_g.to_csv(top100, sep = "\t")
	
	top1000_g = counts.iloc[:1000, :]
	top1000_g.to_csv(top1000, sep = "\t")
	
	# Generate chaperones and variants dataframe
	# Generate chaperones and non replicative variants dataframe
	Tools.create_folder(("/").join(str(full).split("/")[:-1]))
	Tools.create_folder(("/").join(str(nonRcv).split("/")[:-1]))

	if Config.args.which == "variants_chaperones":
		for i in counts.index.to_list():
			counts.rename(index = {i: i.split('.')[0]}, inplace = True)
		
		cv_df = counts.loc[chaperones_variants.index.to_list()]
		nonRv_df = counts.loc[chaperone_nonRv.index.to_list()]	

		cv_df.to_csv(str(full), sep = "\t")
		nonRv_df.to_csv(str(nonRcv), sep = "\t")
	
	# Generate random sets
	Tools.create_folder(rand)

	for c in range(10):
		rng = default_rng()
		r = rng.choice(len(counts), size = 125, replace = False)
		df_random = pd.DataFrame(columns = counts.columns)
		
		for i in r:
			df_random = df_random.append(counts.iloc[i])

		df_random.to_csv(rand.joinpath(str(c) + ".tsv"), sep = "\t")
	
	# Generate counts by tissue/cancer
	Tools.create_folder(tissue_counts)

	if Config.args.which == "variants_chaperones":
		counts_by_tissue = cv_df.T
	elif Config.args.which == "Normalized":
		counts_by_tissue = counts.T
	counts_by_tissue = counts_by_tissue.join(metadata[Config.args.smts])
	tissue_types = pd.unique(metadata[Config.args.smts])

	for t in tissue_types:
		df = pd.DataFrame(columns = counts_by_tissue.columns)
		df = df.append(counts_by_tissue[counts_by_tissue[Config.args.smts] == t])
		df.pop(Config.args.smts)
		df = df.T
		df.to_csv(tissue_counts.joinpath(t + ".tsv"), sep = '\t')

	if Config.args.dataset == "GTEx":
		Tools.create_folder(("/").join(str(normal).split("/")[:-1]))
		Tools.create_folder(("/").join(str(wo_bbbpst_tissues).split("/")[:-1]))

		# Generate counts without transformed cells
		if Config.args.which == "variants_chaperones":
			counts_wo_tcells = cv_df.T
		elif Config.args.which == "Normalized":
			counts_wo_tcells = counts.T
		counts_wo_tcells = counts_wo_tcells.join(metadata[Config.args.smts])
		counts_wo_tcells = counts_wo_tcells.join(metadata[Config.args.smtsd])
		tissue_types = list(pd.unique(metadata[Config.args.smtsd]))

		toDiscard = ["Cells - EBV-transformed lymphocytes",
					"Cells - Leukemia cell line (CML)",
					"Cells - Transformed fibroblasts"]

		for discard in toDiscard:
			tissue_types.remove(discard)

		df = pd.DataFrame(columns = counts_wo_tcells.columns)
		for t in tissue_types:
			df = df.append(counts_wo_tcells[counts_wo_tcells[Config.args.smtsd] == t])
			
		df.pop(Config.args.smts)
		df.pop(Config.args.smtsd)
		df = df.T
		print(df.head())
		df.to_csv(str(normal), sep = '\t')

		# Generate counts without Brain, Blood, Bone Marrow, Pituitary, Spleen, Testis
		if Config.args.which == "variants_chaperones":
			counts_wo_bbbpst = cv_df.T
		elif Config.args.which == "Normalized":
			counts_wo_bbbpst = counts.T
		counts_wo_bbbpst = counts_wo_bbbpst.join(metadata[Config.args.smts])
		tissue_types = list(pd.unique(metadata[Config.args.smts]))

		toDiscard = ["Brain", "Blood", "Bone Marrow", 
			"Pituitary", "Spleen", "Testis"]

		for discard in toDiscard:
			tissue_types.remove(discard)

		df = pd.DataFrame(columns = counts_wo_bbbpst.columns)
		for t in tissue_types:
			df = df.append(counts_wo_bbbpst[counts_wo_bbbpst[Config.args.smts] == t])
			
		df.pop(Config.args.smts)
		df = df.T
		df.to_csv(str(wo_bbbpst_tissues), sep = '\t')

if __name__ == '__main__':

	generate_data(
		full = Config.args.full,
		meta = Config.args.meta,
		rand = Config.args.rand,
		top100 = Config.args.top100,
		cv_list = Config.args.list,
		nonRcv = Config.args.nonRcv,
		top1000 = Config.args.top1000,
		normal = Config.args.onlyNormal,
		tissue_counts = Config.args.tissue,
		normalized_counts = Config.args.norm,
		nonRcv_list = Config.args.nonReplicative,
		wo_bbbpst_tissues = Config.args.WoTissues)