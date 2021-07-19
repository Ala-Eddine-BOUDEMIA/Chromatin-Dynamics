import pandas as pd 
import plotly.express as px

from pathlib import Path

import Tools
import Config

def explore_data(
	meta, counts, rand,
	qc_imgs, qc_htmls,
	qc_rand_img, qc_rand_html):
	
	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)
	
	for i in range(len(rand_files)):
		link_img = qc_rand_img.joinpath("random" + str(i))
		link_p = qc_rand_html.joinpath("random" + str(i))
		
		Tools.create_folder(link_img)
		Tools.create_folder(link_p)
		
		qc_imgs.append(link_img)
		qc_htmls.append(link_p)

	# read metadata file
	metadata = pd.read_csv(meta, header = 0, sep = "\t")

	# Number of samples per tissue 
	# smtsd or  gdc_cases.tissue_source_site.project
	samples_per_tissue = metadata.groupby(["smtsd"]).agg({'run':'count'})
	samples_per_tissue.sort_values("smtsd")

	# Runs per tissue
	runs_per_tissue = metadata.groupby(["smtsd"]).agg({'run':'unique'})
	runs_per_tissue.sort_values("smtsd")

	# Process each file
	for file, image, html in zip(counts, qc_imgs, qc_htmls):
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
		counts = Config.counts,
		rand = Config.args.rand,
		qc_imgs = Config.general_qc_imgs,
		qc_htmls = Config.general_qc_htmls,
		qc_rand_img = Config.args.IgeneralRand,
		qc_rand_html = Config.args.PgeneralRand)