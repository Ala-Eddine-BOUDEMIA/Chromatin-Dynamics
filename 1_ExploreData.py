import pandas as pd 
import plotly.express as px

from pathlib import Path

import Config

def explore_data(
	meta, #raw_counts, filtered_counts, 
	counts_norm, top1000, 
	cv_counts, #rand, 
	generalNormImages, #generalRandImages,
	generalTop1000Images, generalCvImages, generalNormPlotly, 
	#generalRandPlotly, 
	generalTop1000Ploty, generalCvPlotly):
	
	counts = [counts_norm, top1000, cv_counts]
	images = [generalNormImages, generalTop1000Images, generalCvImages]
	htmls = [generalNormPlotly, generalTop1000Ploty, generalCvPlotly]

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
			fig.write_image(str(image.joinpath(t.replace(" ", "_") + ".png")))
			fig.show()

if __name__ == '__main__':	
	
	explore_data(
		meta = Config.args.meta,
		#raw_counts = Config.args.bf,
		#filtered_counts = Config.args.af,
		counts_norm = Config.args.norm,
		top1000 = Config.args.top1000,
		cv_counts = Config.args.cv,
		#rand = Config.args.rand,
		generalNormImages = Config.args.IgeneralNorm,
		#generalRandImages = Config.args.IgeneralRand,
		generalTop1000Images = Config.args.IgeneralTop,
		generalCvImages = Config.args.IgeneralCV,
		generalNormPlotly = Config.args.PgeneralNorm,
		#generalRandPlotly = Config.args.PgeneralRand,
		generalTop1000Ploty = Config.args.PgeneralTop,
		generalCvPlotly = Config.args.PgeneralCV)