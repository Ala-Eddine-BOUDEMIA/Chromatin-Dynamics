import numpy as np
import pandas as pd
import plotly.express as px

import Config

def mean_std(
	raw_counts, filtered_counts, counts_norm, top1000, cv_counts, rand, raw_mv_i, 
	filtered_mv_i, counts_norm_mv_image, rand_images, top1000_mv_image, cv_mv_image, 
	counts_norm_mv_p, rand_p, top1000_mv_p, cv_mv_p, raw_mv_p, filtered_mv_p):
	
	counts = [raw_counts, filtered_counts, counts_norm, top1000, cv_counts]
	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)

	images = [raw_mv_i, filtered_mv_i, counts_norm_mv_image, top1000_mv_image, cv_mv_image]
	htmls = [raw_mv_p, filtered_mv_p, counts_norm_mv_p, top1000_mv_p, cv_mv_p]

	for i in range(len(rand_files)):
		images.append(rand_images.joinpath("random" + str(i) + ".png"))
		htmls.append(rand_p.joinpath("random" + str(i) + ".html"))

	for file, image, html in zip(counts, images, htmls):
		f = pd.read_csv(file, header = 0, index_col = 0, sep = "\t")

		mean = np.log2(f + 0.5).mean(axis = 1)
		std = np.sqrt(f).std(axis = 1)
		meanStd = pd.DataFrame(data=[mean, std], 
			index=['log2(mean+0.5)','sqrt(std)']).T

		fig = px.scatter(
			data_frame = meanStd, 
			x = "log2(mean+0.5)", 
			y = "sqrt(std)",
			hover_data = [meanStd.index],
			title = "Mean Variance Plot")

		fig.write_html(str(html))
		fig.write_image(str(image), width = 2048, height = 1024)
		#fig.show()

if __name__ == '__main__':
	
	mean_std(
		raw_counts = Config.args.bf,
		filtered_counts = Config.args.af,
		counts_norm = Config.args.norm,
		top1000 = Config.args.top1000,
		cv_counts = Config.args.cv,
		rand = Config.args.rand,
		raw_mv_i = Config.args.ImvRaw,
		filtered_mv_i = Config.args.ImvFiltered,
		counts_norm_mv_image = Config.args.ImvNorm,
		rand_images = Config.args.ImvRand,
		top1000_mv_image = Config.args.ImvTop,
		cv_mv_image = Config.args.ImvCV,
		counts_norm_mv_p = Config.args.PmvNorm,
		rand_p = Config.args.PmvRand,
		top1000_mv_p = Config.args.PmvTop,
		cv_mv_p = Config.args.PmvCV,
		raw_mv_p = Config.args.PmvRaw,
		filtered_mv_p = Config.args.PmvFiltered)