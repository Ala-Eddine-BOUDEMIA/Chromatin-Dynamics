import gc
import numpy as np
import pandas as pd

import Tools
import Config 

def correlation(
	counts, rand, g_corr, s_corr,
	by_tissues, g_corr_rand, s_corr_rand,
	g_corr_by_tissue, s_corr_by_tissue):

	for i, j in zip(g_corr, s_corr):
		Tools.create_folder(('/').join(str(i)).split("/")[:-1])
        Tools.create_folder(('/').join(str(j)).split("/")[:-1])
	
	tissue_files = sorted([f for f in by_tissues.glob("**/*.tsv") if f.is_file()])
	for path in tissue_files:
		counts.append(path) 
	
	rand_files = sorted([f for f in rand.glob("**/*.tsv") if f.is_file()])
	for path in rand_files:
		counts.append(path)

	for i in range(len(tissue_files)):
		tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
		link_g = g_corr_by_tissue.joinpath(tissue_name)
		link_s = s_corr_by_tissue.joinpath(tissue_name)
		
		Tools.create_folder(link_g)
		Tools.create_folder(link_s)
		
		g_corr.append(link_g.joinpath(tissue_name + ".tsv"))
		s_corr.append(link_s.joinpath(tissue_name + ".tsv"))
		
	Tools.create_folder(g_corr_rand)
	Tools.create_folder(s_corr_rand)
	for i in range(len(rand_files)):
		g_corr.append(g_corr_rand.joinpath("random" + str(i) + ".tsv"))
		s_corr.append(s_corr_rand.joinpath("random" + str(i) + ".tsv"))

	# Process each file
	for filee, g, s in zip(counts, g_corr, s_corr):
		f = pd.read_csv(filee, header = 0, index_col = 0, sep = "\t")
		f_log2 = pd.DataFrame(np.log2(f + 1))

		# Samples
		s_f_log2 = f_log2.corr()
		s_f_log2.to_csv(s, sep = "\t", float_format='%.3f')

		# Genes
		g_f_log2 = f_log2.T.corr()
		g_f_log2.to_csv(g, sep = "\t", float_format='%.3f')

		# Free memory
		del(f)
		del(g_f_log2)
		del(s_f_log2)
		gc.collect()

if __name__ == '__main__':

	correlation(
		counts = Config.counts,
		rand = Config.args.rand,
		g_corr = Config.g_corr,
		s_corr = Config.s_corr,
		by_tissues = Config.args.tissue,
		g_corr_rand = Config.args.corrRandG,
		s_corr_rand = Config.args.corrRandS,
		g_corr_by_tissue = Config.args.corrTissueG,
		s_corr_by_tissue = Config.args.corrTissueS)