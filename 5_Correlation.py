import numpy as np
import pandas as pd

import Tools
import Config 

def correlation(
	rand, top88, cv_counts, top1000, nrcv_counts,
	g_corr_cv , s_corr_cv, g_corr_top88, s_corr_top88,
	g_corr_rand, s_corr_rand, g_corr_top1000, s_corr_top1000, 
	g_corr_nrcv, s_corr_nrcv):
	
	counts = [top1000, top88, cv_counts, nrcv_counts] 

	
	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)
	
	g_corr = [] #[g_corr_top1000, g_corr_top88, g_corr_cv, g_corr_nrcv]

	s_corr = [] #[s_corr_top1000, s_corr_top88, s_corr_cv, s_corr_nrcv] 
	
	for i in range(len(rand_files)):
		g_corr.append(g_corr_rand.joinpath("random" + str(i) + ".tsv"))
		s_corr.append(s_corr_rand.joinpath("random" + str(i) + ".tsv"))

	# Process each file
	for file, g, s in zip(counts, g_corr, s_corr):
		f = pd.read_csv(file, header = 0, index_col = 0, sep = "\t")
		f_log2 = pd.DataFrame(np.log2(f + 1))

		# Samples
		s_f_log2 = f_log2.corr('pearson')
		s_f_log2.to_csv(s, sep = "\t")

		# Genes
		g_f_log2 = f_log2.T.corr('pearson')
		g_f_log2.to_csv(g, sep = "\t")

if __name__ == '__main__':

	correlation(
		rand = Config.args.rand,
		top88 = Config.args.top88,
		cv_counts = Config.args.cv,
		top1000 = Config.args.top1000,
		nrcv_counts = Config.args.nonRcv,
		g_corr_cv = Config.args.corrCVg,
		s_corr_cv = Config.args.corrCVs,
		g_corr_top88 = Config.args.corr88G,
		s_corr_top88 = Config.args.corr88S,
		g_corr_rand = Config.args.corrRandG,
		s_corr_rand = Config.args.corrRandS,
		g_corr_top1000 = Config.args.corrTopG,
		s_corr_top1000 = Config.args.corrTopS,
		g_corr_nrcv = Config.args.corrNonRcvG,
		s_corr_nrcv = Config.args.corrNonRcvS)