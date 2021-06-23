import pandas as pd 
import plotly.express as px

from pathlib import Path

import Tools
import Config

def explore_data(
	meta, raw_counts, filtered_counts, counts_norm, top1000, cv_counts, rand, 
	generalRawImages, generalFilteredImages, generalNormImages, generalRandImages,
	generalTop1000Images, generalCvImages, generalNormPlotly, generalRandPlotly, 
	generalTop1000Plotly, generalCvPlotly, generalRawPlotly, generalFilteredPlotly):
	
	counts = [raw_counts, filtered_counts, counts_norm, top1000, cv_counts]
	
	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)

	images = [generalRawImages, generalFilteredImages, generalNormImages, generalTop1000Images, generalCvImages]
	htmls = [generalRawPlotly, generalFilteredPlotly, generalNormPlotly, generalTop1000Plotly, generalCvPlotly]
	
	for i in range(len(rand_files)):
		link_img = generalRandImages.joinpath("random" + str(i))
		link_p = generalRandPlotly.joinpath("random" + str(i))
		
		Tools.create_folder(link_img)
		Tools.create_folder(link_p)
		
		images.append(link_img)
		htmls.append(link_p)

	# read metadata file
	metadata = pd.read_csv(meta, header = 0, sep = "\t")

	# Number of samples per tissue
	samples_per_tissue = metadata.groupby(["smtsd"]).agg({'run':'count'})
	samples_per_tissue.sort_values("smtsd")

	# Runs per tissue
	runs_per_tissue = metadata.groupby(["smtsd"]).agg({'run':'unique'})
	runs_per_tissue.sort_values("smtsd")

	# Process each file
	for file, image, html in zip(counts, images, htmls):
		f = pd.read_csv(file, header = 0, index_col = 0, sep = "\t")

		# Number of counts per sample
		counts_per_sample = f.sum(axis = 0)

		# Number of counts per gene
		counts_per_gene = f.sum(axis = 1)

		# cpm
		cpm = f
		total = counts_per_sample.div(1e6)
		cpm = cpm.loc[:,:].div(total) 
		cpm_sum = cpm.sum(axis=1)

		# Number of counts per tissue
		counts_per_tissue_sample = []		
		counts_per_tissue_gene = []

		for tissue in runs_per_tissue.index:
			for runs in runs_per_tissue.loc[tissue]:
				# per sample
				counts_per_sample_per_tissue = f.loc[:,runs].sum(axis=0)
				counts_per_tissue_sample.append([list(counts_per_sample_per_tissue), tissue])

				# per gene
				counts_per_gene_per_tissue = f.loc[:,runs].sum(axis=1)
				counts_per_tissue_gene.append((counts_per_gene_per_tissue, tissue))
				
		# per sample
		df_counts_per_tissue_sample = pd.DataFrame(
			[x[0] for x in counts_per_tissue_sample], 
			index = [x[1] for x in counts_per_tissue_sample])

		# per gene
		df_counts_per_tissue_gene = pd.DataFrame(
			[x[0] for x in counts_per_tissue_gene], 
			index = [x[1] for x in counts_per_tissue_gene])

		dataframes = [
			pd.DataFrame(counts_per_sample, columns = ["Counts Per Sample"]), 
			pd.DataFrame(counts_per_gene, columns = ["Counts Per Gene"]), 
			pd.DataFrame(cpm_sum, columns = ["Counts Per Million"]), 
			df_counts_per_tissue_sample.T, 
			df_counts_per_tissue_gene.T]

		titles = [
			"Counts per sample", "Counts per gene", "Counts per million",
			"Counts per sample by tissue", "Counts per gene by tissue"]

		for d, t in zip(dataframes, titles):
			fig = px.box(
				data_frame = d,
				y = d.columns,
				title = t)

			fig.write_html(str(html.joinpath(t.replace(" ", "_") + ".html")))
			fig.write_image(str(image.joinpath(t.replace(" ", "_") + ".png")), 
				width = 2048, height = 1024)
			fig.show()

if __name__ == '__main__':	
	
	explore_data(
		meta = Config.args.meta,
		raw_counts = Config.args.bf,
		filtered_counts = Config.args.af,
		counts_norm = Config.args.norm,
		top1000 = Config.args.top1000,
		cv_counts = Config.args.cv,
		rand = Config.args.rand,
		generalRawImages = Config.args.IgeneralRaw,
		generalFilteredImages = Config.args.IgeneralFiltered,
		generalNormImages = Config.args.IgeneralNorm,
		generalRandImages = Config.args.IgeneralRand,
		generalTop1000Images = Config.args.IgeneralTop,
		generalCvImages = Config.args.IgeneralCV,
		generalNormPlotly = Config.args.PgeneralNorm,
		generalRandPlotly = Config.args.PgeneralRand,
		generalTop1000Plotly = Config.args.PgeneralTop,
		generalCvPlotly = Config.args.PgeneralCV,
		generalRawPlotly = Config.args.PgeneralRaw,
		generalFilteredPlotly = Config.args.PgeneralFiltered)