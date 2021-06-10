import Plots
import pandas as pd 

gtex_counts = pd.read_csv(
	"Data/GTEx/read_counts.tsv", 
	header=0, index_col=0, sep="\t")

gtex = pd.read_csv(
	"Data/GTEx/GTEx.tsv",
	header=0, sep="\t")

# 1- Number of samples per tissue
samples_per_tissue = gtex.groupby(["smts"]).agg({'run':'count'})

Plots.histogram(samples_per_tissue, 
				title="GTEx Number of Samples Per Tissue",
				x=samples_per_tissue.index, 
				y=samples_per_tissue.columns)

Plots.boxplot(samples_per_tissue, 
				title="GTEx Number of Samples Per Tissue",
				x=samples_per_tissue.columns)

# 2- Number of counts per sample
counts_per_sample = gtex_counts.sum(axis=0)

Plots.histogram(counts_per_sample, 
				title="GTEx Number of Counts Per Sample",
				x=counts_per_sample)

Plots.boxplot(counts_per_sample, 
				title="GTEx Number of Counts Per Sample",
				x=counts_per_sample)

# 3- Mean count per sample
mean_counts_per_sample = gtex_counts.mean(axis=0)

Plots.histogram(mean_counts_per_sample, 
				title="GTEx Mean Counts Per Sample",
				x=mean_counts_per_sample)

Plots.boxplot(mean_counts_per_sample, 
				title="GTEx Mean Counts Per Sample",
				x=mean_counts_per_sample)

# 4- Number of counts per gene
counts_per_gene = gtex_counts.sum(axis=1)

Plots.boxplot(counts_per_gene, 
				title="GTEx Number of Counts Per Gene",
				x=counts_per_gene)

# 5- Mean count per gene
mean_counts_per_gene = gtex_counts.mean(axis=1)

Plots.boxplot(mean_counts_per_gene, 
				title="GTEx Mean Counts Per Gene",
				x=mean_counts_per_gene)

# 6- Number of counts per tissue
runs_per_tissue = gtex.groupby(["smts"]).agg({'run':'unique'})

counts_per_tissue_sample, mean_counts_per_tissue_sample = [], []
counts_per_tissue_gene, mean_counts_per_tissue_gene = [], []

for tissue in runs_per_tissue.index:
	for runs in runs_per_tissue.loc[tissue]:
		try:
			# per sample
			counts_per_sample_per_tissue = gtex_counts.loc[:,runs].sum(axis=0)
			counts_per_tissue_sample.append((counts_per_sample_per_tissue, tissue))

			# mean per sample
			mean_counts_per_sample_per_tissue = gtex_counts.loc[:,runs].mean(axis=0)
			mean_counts_per_tissue_sample.append((mean_counts_per_sample_per_tissue, tissue))

			# per gene
			counts_per_gene_per_tissue = gtex_counts.loc[:,runs].sum(axis=1)
			counts_per_tissue_gene.append((counts_per_gene_per_tissue, tissue))
			# mean per gene
			mean_counts_per_gene_per_tissue = gtex_counts.loc[:,runs].mean(axis=1)
			mean_counts_per_tissue_gene.append((mean_counts_per_gene_per_tissue, tissue))
		except:
			pass

Plots.histogram(counts_per_tissue_sample, 
				title="GTEx counts_per_tissue_sample",
				x=[x[1] for x in counts_per_tissue_sample],
				y=[y[0].sum() for y in counts_per_tissue_sample])