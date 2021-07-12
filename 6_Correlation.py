import gc
import numpy as np
import pandas as pd

import Tools
import Config 

def correlation(
	rand, top88, cv_counts, top1000, nrcv_counts, normal, wotissues, 
	by_tissues, g_corr_cv , s_corr_cv, g_corr_top88, s_corr_top88,
	g_corr_rand, s_corr_rand, g_corr_top1000, s_corr_top1000, 
	g_corr_nrcv, s_corr_nrcv, g_corr_normal, s_corr_normal,
	g_corr_by_tissue, s_corr_by_tissue, g_corr_wo_tissues, 
	s_corr_wo_tissues):
	
	counts = [normal, wotissues, top1000, top88, cv_counts, nrcv_counts]

	tissue_files = sorted([f for f in by_tissues.iterdir() if f.is_file()])
	for path in tissue_files:
		counts.append(path) 

	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)
	
	g_corr = [g_corr_normal, g_corr_wo_tissues, 
		g_corr_top1000, g_corr_top88, g_corr_cv, g_corr_nrcv]

	s_corr = [s_corr_normal, s_corr_wo_tissues,
		s_corr_top1000, s_corr_top88, s_corr_cv, s_corr_nrcv] 

	for i in range(len(tissue_files)):
		tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
		link_g = g_corr_by_tissue.joinpath(tissue_name)
		link_s = s_corr_by_tissue.joinpath(tissue_name)
		
		Tools.create_folder(link_g)
		Tools.create_folder(link_s)
		
		g_corr.append(link_g.joinpath(tissue_name + ".png"))
		s_corr.append(link_s.joinpath(tissue_name + ".html"))
	
	for i in range(len(rand_files)):
		g_corr.append(g_corr_rand.joinpath("random" + str(i) + ".tsv"))
		s_corr.append(s_corr_rand.joinpath("random" + str(i) + ".tsv"))

	# Process each file
	for filee, g, s in zip(counts, g_corr, s_corr):
		f = pd.read_csv(filee, header = 0, index_col = 0, sep = "\t")
		f_log2 = pd.DataFrame(np.log2(f + 1))

		# Samples
		s_f_log2 = f_log2.corr('pearson')
		s_f_log2.to_csv(s, sep = "\t")

		# Genes
		g_f_log2 = f_log2.T.corr('pearson')
		g_f_log2.to_csv(g, sep = "\t")

		# Free memory
		del(f)
		del(g_f_log2)
		del(s_f_log2)
		gc.collect()

if __name__ == '__main__':

	correlation(
		rand = Config.args.rand,
		top88 = Config.args.top88,
		cv_counts = Config.args.cv,
		top1000 = Config.args.top1000,
		nrcv_counts = Config.args.nonRcv,
		normal = Config.args.onlyNormal,
		wotissues = Config.args.WoTissues,
		by_tissues = Config.args.tissue,
		g_corr_cv = Config.args.corrCVg,
		s_corr_cv = Config.args.corrCVs,
		g_corr_top88 = Config.args.corr88G,
		s_corr_top88 = Config.args.corr88S,
		g_corr_rand = Config.args.corrRandG,
		s_corr_rand = Config.args.corrRandS,
		g_corr_top1000 = Config.args.corrTopG,
		s_corr_top1000 = Config.args.corrTopS,
		g_corr_nrcv = Config.args.corrNonRcvG,
		s_corr_nrcv = Config.args.corrNonRcvS,
		g_corr_normal = Config.args.corrNormalG,
		s_corr_normal = Config.args.corrNormalS,
		g_corr_by_tissue = Config.args.corrTissueG,
		s_corr_by_tissue = Config.args.corrTissueS,
		g_corr_wo_tissues = Config.args.corrWoTissuesG,
		s_corr_wo_tissues = Config.args.corrWoTissuesS)