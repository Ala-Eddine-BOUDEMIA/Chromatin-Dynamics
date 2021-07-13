import gc
import numpy as np
import pandas as pd
import plotly.express as px

import Tools
import Config

def mean_std(
	rand, top100, cv_counts, raw_counts, cv_mv_p, rand_p, 
	top1000, raw_mv_i, raw_mv_p, counts_norm, normal,	
	cv_mv_image, filtered_counts, top1000_mv_p, top100_mv_p,
	rand_images, tissue_counts, normal_mv_p, tissues_mv_p, 
	normal_mv_img, top100_mv_image, top1000_mv_image, 
	counts_norm_mv_p, tissues_images, filtered_mv_p, 
	filtered_mv_i, counts_norm_mv_image, counts_without_tissues,
	counts_without_tissues_mv_p, counts_without_tissues_mv_image):
	
	counts = [counts_norm, top1000, top100, cv_counts] 
	#raw_counts, filtered_counts, normal, counts_without_tissues,
	
	tissue_files = sorted([f for f in tissue_counts.iterdir() if f.is_file()])
	for path in tissue_files:
		counts.append(path)

	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)

	images = [counts_norm_mv_image, top1000_mv_image, top100_mv_image, cv_mv_image]
	#raw_mv_i, filtered_mv_i,  normal_mv_img, counts_without_tissues_mv_image, 

	htmls = [counts_norm_mv_p, top1000_mv_p, top100_mv_p, cv_mv_p]
	#raw_mv_p, filtered_mv_p, normal_mv_p, counts_without_tissues_mv_p, 
	
	for i in range(len(tissue_files)):
		tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
		link_img = tissues_images.joinpath(tissue_name)
		link_p = tissues_mv_p.joinpath(tissue_name)

		Tools.create_folder(link_img)
		Tools.create_folder(link_p)

		images.append(link_img.joinpath(tissue_name + ".png"))
		htmls.append(link_p.joinpath(tissue_name + ".html"))
	
	for i in range(len(rand_files)):
		images.append(rand_images.joinpath("random" + str(i) + ".png"))
		htmls.append(rand_p.joinpath("random" + str(i) + ".html"))

	for file, image, html in zip(counts, images, htmls):
		f = pd.read_csv(file, header = 0, index_col = 0, sep = "\t")

		mean = np.log2(f + 1).mean(axis = 1)
		std = np.sqrt(f).std(axis = 1)
		meanStd = pd.DataFrame(data=[mean, std], 
			index=['log2(mean+1)','sqrt(std)']).T

		fig = px.scatter(
			data_frame = meanStd, 
			x = "log2(mean+1)", 
			y = "sqrt(std)",
			hover_data = [meanStd.index],
			title = "Mean Variance Plot")

		fig.write_html(str(html))
		fig.write_image(str(image), width = 2048, height = 1024)
		#fig.show()

		del(f)
		del(fig)
		gc.collect()

if __name__ == '__main__':
	
	mean_std(
		rand = Config.args.rand,
		top100 = Config.args.top100,
		cv_counts = Config.args.cv,
		raw_counts = Config.args.bf,
		cv_mv_p = Config.args.PmvCV,
		rand_p = Config.args.PmvRand,
		top1000 = Config.args.top1000,
		raw_mv_i = Config.args.ImvRaw,
		raw_mv_p = Config.args.PmvRaw,
		counts_norm = Config.args.norm,	
		normal = Config.args.onlyNormal,
		cv_mv_image = Config.args.ImvCV,
		filtered_counts = Config.args.af,
		top1000_mv_p = Config.args.PmvTop,
		top100_mv_p = Config.args.PmvTop100,
		rand_images = Config.args.ImvRand,
		tissue_counts = Config.args.tissue,
		normal_mv_p = Config.args.PmvNormal,
		tissues_mv_p = Config.args.PmvTissue,
		normal_mv_img = Config.args.ImvNormal,
		top100_mv_image = Config.args.ImvTop100,
		top1000_mv_image = Config.args.ImvTop,
		counts_norm_mv_p = Config.args.PmvNorm,
		tissues_images = Config.args.ImvTissue,
		filtered_mv_p = Config.args.PmvFiltered,
		filtered_mv_i = Config.args.ImvFiltered,
		counts_norm_mv_image = Config.args.ImvNorm,
		counts_without_tissues = Config.args.WoTissues,
		counts_without_tissues_mv_p = Config.args.PmvWoTissues,
		counts_without_tissues_mv_image = Config.args.ImvWoTissues)