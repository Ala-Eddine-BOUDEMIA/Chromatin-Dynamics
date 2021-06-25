import numpy as np
import pandas as pd

import Config 

def correlation(
	counts_norm, top1000, top76, 
	cv_counts, nrcv_counts, rand, 
	g_corr_norm, g_corr_top1000, g_corr_top76, 
	g_corr_cv, g_corr_nrcv, g_corr_rand,
	s_corr_norm, s_corr_top1000, s_corr_top76, 
	s_corr_rand, s_corr_nrcv, s_corr_cv):
	
	counts = [counts_norm, top1000, top76, cv_counts, nrcv_counts] 
	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)

	g_corr = [g_corr_norm, g_corr_top1000, g_corr_top76, g_corr_cv, g_corr_nrcv] 
	s_corr = [s_corr_norm, s_corr_top1000, s_corr_top76, s_corr_cv, s_corr_nrcv]

	for i in range(len(rand_files)):
		g_corr.append(g_corr_rand.joinpath("random" + str(i) + ".tsv"))
		s_corr.append(s_corr_rand.joinpath("random" + str(i) + ".tsv"))

	# Process each file
	for file, g, s in zip(counts, g_corr, s_corr):
		f = pd.read_csv(file, header = 0, index_col = 0, sep = "\t")
		f_log2 = pd.DataFrame(np.log2(f + 1))

		# Samples
		s_f_log2 = f_log2.corr()
		s_f_log2.to_csv(s, sep = "\t")

		# Genes
		g_f_log2 = f_log2.T.corr()
		g_f_log2.to_csv(g, sep = "\t")

if __name__ == '__main__':

	correlation(
		counts_norm = Config.args.norm,
		top1000 = Config.args.top1000,
		top76 = Config.args.top76,
		cv_counts = Config.args.cv,
		nrcv_counts = Config.args.nonRcv,
		rand = Config.args.rand,
		g_corr_norm = Config.args.corrNormG,
		g_corr_top1000 = Config.args.corrTopG,
		g_corr_top76 = Config.args.corr76G,
		g_corr_cv = Config.args.corrCVg,
		g_corr_nrcv = Config.args.corrNonRcvG,
		g_corr_rand = Config.args.corrRandG,
		s_corr_norm = Config.args.corrNormS,
		s_corr_top1000 = Config.args.corrTopS,
		s_corr_top76 = Config.args.corr76S,
		s_corr_rand = Config.args.corrRandS,
		s_corr_cv = Config.args.corrCVs,
		s_corr_nrcv = Config.args.corrNonRcvS)