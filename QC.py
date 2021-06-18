import Plots

import pandas as pd 

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

gtex_counts_before_filtering = pd.read_csv(
	"Data/GTEx/BeforeFiltering/PairedEndRounded.tsv", 
	header=0, index_col=0, sep="\t")

gtex_counts_after_filtering = pd.read_csv(
	"Data/GTEx/AfterFiltering/GTExFilteredCPM10S36.tsv", 
	header=0, index_col=0, sep="\t")

gtex = pd.read_csv(
	"Data/GTEx/Metadata/GTEx.tsv",
	header=0, sep="\t")

# 1- Number of samples per tissue
samples_per_tissue = gtex.groupby(["smtsd"]).agg({'run':'count'})
samples_per_tissue.sort_values("smtsd")

Plots.histogram(samples_per_tissue, 
				title="GTEx Number of Samples Per Tissue",
				x=samples_per_tissue.index, 
				y=samples_per_tissue.columns)

Plots.boxplot(samples_per_tissue, 
				title="GTEx Number of Samples Per Tissue",
				x=samples_per_tissue.columns)

# 2- Number of counts per sample
# 2.1- Before Filtering
counts_per_sample_before_filtering = gtex_counts_before_filtering.sum(axis=0)

Plots.boxplot(counts_per_sample_before_filtering, 
				title="GTEx Number of Counts Per Sample Before Filtering",
				x=counts_per_sample_before_filtering)

# 2.2- After Filtering
counts_per_sample_after_filtering = gtex_counts_after_filtering.sum(axis=0)

Plots.boxplot(counts_per_sample_after_filtering, 
				title="GTEx Number of Counts Per Sample After Filtering",
				x=counts_per_sample_after_filtering)

# 3- Mean count per sample
# 3.1- Before Filtering
mean_counts_per_sample_before_filtering = gtex_counts_before_filtering.mean(axis=0)

Plots.boxplot(mean_counts_per_sample_before_filtering, 
				title="GTEx Mean Counts Per Sample Before Filtering",
				x=mean_counts_per_sample_before_filtering)

# 3.1- After Filtering
mean_counts_per_sample_after_filtering = gtex_counts_after_filtering.mean(axis=0)

Plots.boxplot(mean_counts_per_sample_after_filtering, 
				title="GTEx Mean Counts Per Sample After Filtering",
				x=mean_counts_per_sample_after_filtering)

# 4- Number of counts per gene
# 4.1- Before Filtering
counts_per_gene_before_filtering = gtex_counts_before_filtering.sum(axis=1)

Plots.boxplot(counts_per_gene_before_filtering , 
				title="GTEx Number of Counts Per Gene Before Filtering",
				x=counts_per_gene_before_filtering )

# 4.2- After Filtering
counts_per_gene_after_filtering = gtex_counts_after_filtering.sum(axis=1)

Plots.boxplot(counts_per_gene_after_filtering , 
				title="GTEx Number of Counts Per Gene After Filtering",
				x=counts_per_gene_after_filtering )

# 5- Mean count per gene
# 5.1- Before Filtering
mean_counts_per_gene_before_filtering = gtex_counts_before_filtering.mean(axis=1)

Plots.boxplot(mean_counts_per_gene_before_filtering, 
				title="GTEx Mean Counts Per Gene Before Filtering",
				x=mean_counts_per_gene_before_filtering)

# 5.2- After Filtering
mean_counts_per_gene_after_filtering = gtex_counts_after_filtering.mean(axis=1)

Plots.boxplot(mean_counts_per_gene_after_filtering, 
				title="GTEx Mean Counts Per Gene After Filtering",
				x=mean_counts_per_gene_after_filtering)

# 6- CPM
# 6.1- CPM Before Filtering
cpm_before_filtering = gtex_counts_before_filtering
total_before_filtering  = counts_per_sample_before_filtering.div(1e6)
cpm_before_filtering  = cpm_before_filtering.loc[:,:].div(total_before_filtering ) 
cpm_before_filtering = cpm_before_filtering.sum(axis=1)

Plots.boxplot(cpm_before_filtering, title="CPM Before Filtering")

# 6.2- CPM After Filtering
cpm_after_filtering = gtex_counts_after_filtering
total_after_filtering = counts_per_sample_after_filtering.div(1e6)
cpm_after_filtering= cpm_after_filtering.loc[:,:].div(total_after_filtering) 
cpm_after_filtering = cpm_after_filtering.sum(axis=1)

Plots.boxplot(cpm_after_filtering, title="CPM After Filtering")

# 7- Number of counts per tissue
runs_per_tissue = gtex.groupby(["smtsd"]).agg({'run':'unique'})
runs_per_tissue.sort_values("smtsd")

# 7.1- Before Filtering
counts_per_tissue_sample_before_filtering, mean_counts_per_tissue_sample_before_filtering = [], []
counts_per_tissue_gene_before_filtering, mean_counts_per_tissue_gene_before_filtering = [], []

for tissue in runs_per_tissue.index:
	for runs in runs_per_tissue.loc[tissue]:
		try:
			# per sample
			counts_per_sample_per_tissue_before_filtering = gtex_counts_before_filtering.loc[:,runs].sum(axis=0)
			counts_per_tissue_sample_before_filtering.append([list(counts_per_sample_per_tissue_before_filtering), tissue])

			# mean per sample
			mean_counts_per_sample_per_tissue_before_filtering = gtex_counts_before_filtering.loc[:,runs].mean(axis=0)
			mean_counts_per_tissue_sample_before_filtering.append((mean_counts_per_sample_per_tissue_before_filtering, tissue))

			# per gene
			counts_per_gene_per_tissue_before_filtering = gtex_counts_before_filtering.loc[:,runs].sum(axis=1)
			counts_per_tissue_gene_before_filtering.append((counts_per_gene_per_tissue_before_filtering, tissue))
			
			# mean per gene
			mean_counts_per_gene_per_tissue_before_filtering = gtex_counts_before_filtering.loc[:,runs].mean(axis=1)
			mean_counts_per_tissue_gene_before_filtering.append((mean_counts_per_gene_per_tissue_before_filtering, tissue))

		except:
			pass

# per sample
df_counts_per_tissue_sample_before_filtering = pd.DataFrame([x[0] for x in counts_per_tissue_sample_before_filtering], 
									index=[x[1] for x in counts_per_tissue_sample_before_filtering])

Plots.boxplot(df_counts_per_tissue_sample_before_filtering.transpose(), 
				title="GTEx Number of Counts Per Sample Per Tissue Before Filtering")

# mean per sample
df_mean_counts_per_tissue_sample_before_filtering = pd.DataFrame([x[0] for x in mean_counts_per_tissue_sample_before_filtering], 
									index=[x[1] for x in mean_counts_per_tissue_sample_before_filtering])

Plots.boxplot(df_mean_counts_per_tissue_sample_before_filtering.transpose(), 
				title="GTEx Mean Counts Per Sample Per Tissue Before Filtering")

# per gene
df_counts_per_tissue_gene_before_filtering = pd.DataFrame([x[0] for x in counts_per_tissue_gene_before_filtering], 
									index=[x[1] for x in counts_per_tissue_gene_before_filtering])

Plots.boxplot(df_counts_per_tissue_gene_before_filtering.transpose(), 
				title="GTEx Number of Counts Per Gene Per Tissue Before Filtering")

# mean per gene
df_mean_counts_per_tissue_gene_before_filtering = pd.DataFrame([x[0] for x in mean_counts_per_tissue_gene_before_filtering], 
									index=[x[1] for x in mean_counts_per_tissue_gene_before_filtering])

Plots.boxplot(df_mean_counts_per_tissue_gene_before_filtering.transpose(), 
				title="GTEx Mean Counts Per Gene Per Tissue Before Filtering")

# 7.1- After Filtering
counts_per_tissue_sample_after_filtering, mean_counts_per_tissue_sample_after_filtering = [], []
counts_per_tissue_gene_after_filtering, mean_counts_per_tissue_gene_after_filtering = [], []

for tissue in runs_per_tissue.index:
	for runs in runs_per_tissue.loc[tissue]:
		try:
			# per sample
			counts_per_sample_per_tissue_after_filtering = gtex_counts_after_filtering.loc[:,runs].sum(axis=0)
			counts_per_tissue_sample_after_filtering.append([list(counts_per_sample_per_tissue_after_filtering), tissue])

			# mean per sample
			mean_counts_per_sample_per_tissue_after_filtering = gtex_counts_after_filtering.loc[:,runs].mean(axis=0)
			mean_counts_per_tissue_sample_after_filtering.append((mean_counts_per_sample_per_tissue_after_filtering, tissue))

			# per gene
			counts_per_gene_per_tissue_after_filtering = gtex_counts_after_filtering.loc[:,runs].sum(axis=1)
			counts_per_tissue_gene_after_filtering.append((counts_per_gene_per_tissue_after_filtering, tissue))
			
			# mean per gene
			mean_counts_per_gene_per_tissue_after_filtering = gtex_counts_after_filtering.loc[:,runs].mean(axis=1)
			mean_counts_per_tissue_gene_after_filtering.append((mean_counts_per_gene_per_tissue_after_filtering, tissue))

		except:
			pass

# per sample
df_counts_per_tissue_sample_after_filtering = pd.DataFrame([x[0] for x in counts_per_tissue_sample_after_filtering], 
									index=[x[1] for x in counts_per_tissue_sample_after_filtering])

Plots.boxplot(df_counts_per_tissue_sample_after_filtering.transpose(), 
				title="GTEx Number of Counts Per Sample Per Tissue After Filtering")

# mean per sample
df_mean_counts_per_tissue_sample_after_filtering = pd.DataFrame([x[0] for x in mean_counts_per_tissue_sample_after_filtering], 
									index=[x[1] for x in mean_counts_per_tissue_sample_after_filtering])

Plots.boxplot(df_mean_counts_per_tissue_sample_after_filtering.transpose(), 
				title="GTEx Mean Counts Per Sample Per Tissue After Filtering")

# per gene
df_counts_per_tissue_gene_after_filtering = pd.DataFrame([x[0] for x in counts_per_tissue_gene_after_filtering], 
									index=[x[1] for x in counts_per_tissue_gene_after_filtering])

Plots.boxplot(df_counts_per_tissue_gene_after_filtering.transpose(), 
				title="GTEx Number of Counts Per Gene Per Tissue After Filtering")

# mean per gene
df_mean_counts_per_tissue_gene_after_filtering = pd.DataFrame([x[0] for x in mean_counts_per_tissue_gene_after_filtering], 
									index=[x[1] for x in mean_counts_per_tissue_gene_after_filtering])

Plots.boxplot(df_mean_counts_per_tissue_gene_after_filtering.transpose(), 
				title="GTEx Mean Counts Per Gene Per Tissue After Filtering")