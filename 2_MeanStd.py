import gc
import numpy as np
import pandas as pd
import plotly.express as px

import Tools
import Config

def mean_std(
	counts, mv_p, mv_img,
	rand, rand_p, rand_images,
	tissue_counts, tissues_mv_p, tissues_images):
	
	for i, j in zip(mv_img, mv_p):
		Tools.create_folder(('/').join(str(i)).split("/")[:-1])
		Tools.create_folder(('/').join(str(j)).split("/")[:-1])
		
	tissue_files = sorted([f for f in tissue_counts.glob("**/*.tsv") if f.is_file()])
	for path in tissue_files:
		counts.append(path)

	rand_files = sorted([f for f in rand.glob("**/*.tsv") if f.is_file()])
	for path in rand_files:
		counts.append(path)
	
	for i in range(len(tissue_files)):
		tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
		link_img = tissues_images.joinpath(tissue_name)
		link_p = tissues_mv_p.joinpath(tissue_name)

		Tools.create_folder(link_img)
		Tools.create_folder(link_p)

		mv_img.append(link_img.joinpath(tissue_name + ".png"))
		mv_p.append(link_p.joinpath(tissue_name + ".html"))

	Tools.create_folder(rand_images)
	Tools.create_folder(rand_p)
	for i in range(len(rand_files)):
		mv_img.append(rand_images.joinpath("random" + str(i) + ".png"))
		mv_p.append(rand_p.joinpath("random" + str(i) + ".html"))

	for file, image, html in zip(counts, mv_img, mv_p):
		f = pd.read_csv(file, header = 0, index_col = 0, sep = "\t")

		mean = np.log2(f + 1).mean(axis = 1)
		std = f.std(axis = 1)
		meanStd = pd.DataFrame(data=[mean, std], 
			index=['log2(mean + 1)','std']).T

		fig = px.scatter(
			data_frame = meanStd, 
			x = "log2(mean + 1)", 
			y = "std",
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
		counts = Config.counts,
		mv_p = Config.mv_htmls,
		mv_img = Config.mv_imgs,
		rand = Config.args.rand,
		rand_p = Config.args.PmvRand,
		rand_images = Config.args.ImvRand,
		tissue_counts = Config.args.tissue,
		tissues_mv_p = Config.args.PmvTissue,
		tissues_images = Config.args.ImvTissue)