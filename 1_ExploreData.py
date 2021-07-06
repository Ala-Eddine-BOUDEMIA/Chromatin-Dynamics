import pandas as pd 
import plotly.express as px

from pathlib import Path

import Tools
import Config

def explore_data(
	meta, rand, top88, cv_counts, raw_counts, 
	top1000, counts_norm, normal, filtered_counts, nrcv_counts,
	generalCvImages, generalCvPlotly, generalRawImages, generalRawPlotly, 
	generalNrCvImages, generalNrCvPlotly, generalNormImages, generalNormPlotly, 
	generalRandImages, generalRandPlotly, generalTop88Images, generalTop88Plotly, 
	generalTop1000Images, generalTop1000Plotly, generalNormalImages, generalNormalPlotly,
	generalFilteredImages, generalFilteredPlotly):
	
	counts = [normal]
	"""
	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)"""

	images = [generalNormalImages]
	
	htmls = [generalNormalPlotly]
	"""
	for i in range(len(rand_files)):
		link_img = generalRandImages.joinpath("random" + str(i))
		link_p = generalRandPlotly.joinpath("random" + str(i))
		
		Tools.create_folder(link_img)
		Tools.create_folder(link_p)
		
		images.append(link_img)
		htmls.append(link_p)"""

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
		cpm_sum = cpm.sum(axis = 1)

		# Number of counts per tissue
		counts_per_tissue_sample = []		
		counts_per_tissue_gene = []

		for tissue in runs_per_tissue.index:
			for runs in runs_per_tissue.loc[tissue]:
				try:
					# per sample
					counts_per_sample_per_tissue = f.loc[:, runs].sum(axis = 0)
					counts_per_tissue_sample.append([list(counts_per_sample_per_tissue), tissue])

					# per gene
					counts_per_gene_per_tissue = f.loc[:, runs].sum(axis = 1)
					counts_per_tissue_gene.append((counts_per_gene_per_tissue, tissue))
				except:
					pass
					
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
			#fig.show()

if __name__ == '__main__':	
	
	explore_data(
		meta = Config.args.meta, 
		rand = Config.args.rand,
		top88 = Config.args.top88, 
		cv_counts = Config.args.cv,
		raw_counts = Config.args.bf, 
		top1000 = Config.args.top1000,
		counts_norm = Config.args.norm, 
		normal = Config.args.onlyNormal,
		filtered_counts = Config.args.af,
		nrcv_counts = Config.args.nonRcv, 
		generalCvImages = Config.args.IgeneralCV, 
		generalCvPlotly = Config.args.PgeneralCV,
		generalRawImages = Config.args.IgeneralRaw,
		generalRawPlotly = Config.args.PgeneralRaw,
		generalNrCvImages = Config.args.IgeneralNrCV,
		generalNrCvPlotly = Config.args.PgeneralNrCV,
		generalNormImages = Config.args.IgeneralNorm,
		generalNormPlotly = Config.args.PgeneralNorm,
		generalRandImages = Config.args.IgeneralRand,
		generalRandPlotly = Config.args.PgeneralRand,
		generalTop88Images = Config.args.IgeneralTop88,
		generalTop88Plotly = Config.args.PgeneralTop88,
		generalTop1000Images = Config.args.IgeneralTop,
		generalTop1000Plotly = Config.args.PgeneralTop,
		generalNormalImages = Config.args.IgeneralNormal,
		generalNormalPlotly = Config.args.PgeneralNormal,
		generalFilteredImages = Config.args.IgeneralFiltered,
		generalFilteredPlotly = Config.args.PgeneralFiltered)