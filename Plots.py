import pandas as pd 
import plotly.express as px

# GTEx

gtex = pd.read_csv(
	"Data/GTEx/GTEx.tsv", 
	header=0, sep="\t")

smtsd = px.histogram(gtex, x="smtsd", 
		title="GTEx Number of samples per tissue")
smtsd.show()

smts = px.histogram(gtex, x="smts", 
		title="GTEx Number of samples per tissue")
smts.show()

gtex_counts = pd.read_csv(
	"Data/GTEx/BeforeNorm/gtex_counts_per_sample.tsv", 
	header=0, sep="\t")

gtex_counts_per_sample_hist = px.histogram(gtex_counts, x="Counts",
							title="GTEx Counts per sample")
gtex_counts_per_sample_hist.show()

gtex_counts_per_sample_box = px.box(gtex_counts, x="Counts",
							title="GTEx Counts per sample")
gtex_counts_per_sample_box.show()

gtex_mean_counts = pd.read_csv(
	"Data/GTEx/BeforeNorm/gtex_mean_counts_per_sample.tsv", 
	header=0, sep="\t")

gtex_mean_counts_per_sample_hist = px.histogram(gtex_mean_counts, x="Counts",
									title="GTEx Mean counts per sample")
gtex_mean_counts_per_sample_hist.show()

gtex_mean_counts_per_sample_box = px.box(gtex_mean_counts, x="Counts",
									title="GTEx Mean counts per sample")
gtex_mean_counts_per_sample_box.show()

gtex_counts_tissue = pd.read_csv(
	"Data/GTEx/BeforeNorm/GTEx_CountsPerTissue.tsv", 
	header=0, sep="\t")

gtex_mean_counts_per_tissue_hist = px.histogram(gtex_counts_tissue, 
									x="Tissue", y='Counts',
									title="GTEx Counts per tissue")
gtex_mean_counts_per_tissue_hist.show()

gtex_counts_per_gene = pd.read_csv(
	"Data/GTEx/BeforeNorm/CountsPerGene.tsv",
	header=0, sep="\t")

gtex_counts_per_gene_box = px.box(gtex_counts_per_gene, x="Counts",
							title="GTEx Counts per gene")
gtex_counts_per_gene_box.show()

gtex_counts_per_gene_hist = px.histogram(gtex_counts_per_gene, x="Counts",
							title="GTEx	Counts per gene")
gtex_counts_per_gene_hist.show()

gtex_mean_counts_per_gene = pd.read_csv(
	"Data/GTEx/BeforeNorm/MeanCountsPerGene.tsv",
	header=0, sep="\t")

gtex_mean_counts_per_gene_box = px.box(gtex_mean_counts_per_gene, x="Counts",
								title="GTEx Mean counts per gene")
gtex_mean_counts_per_gene_box.show()

gtex_mean_counts_per_gene_hist = px.histogram(gtex_mean_counts_per_gene, x="Counts",
									title="GTEx Mean counts per gene")
gtex_mean_counts_per_gene_hist.show()

gtex_counts_per_sample_filtered= pd.read_csv(
	"Data/GTEx/BeforeNorm/GTExFilteredCountsPerSample.tsv",
	header=0, sep="\t")

gtex_counts_per_sample_filtered_box = px.box(gtex_counts_per_sample_filtered, x="x",
								title="GTEx counts per sample after filtering")
gtex_counts_per_sample_filtered_box.show()

gtex_counts_per_sample_filtered_hist = px.histogram(gtex_counts_per_sample_filtered, x="x",
									title="GTEx counts per sample after filtering")
gtex_counts_per_sample_filtered_hist.show()

gtex_counts_per_gene_filtered= pd.read_csv(
	"Data/GTEx/BeforeNorm/GTExFilteredCountsPerGene.tsv",
	header=0, sep="\t")

gtex_counts_per_gene_filtered_box = px.box(gtex_counts_per_gene_filtered, x="x",
								title="GTEx counts per gene after filtering")
gtex_counts_per_gene_filtered_box.show()

gtex_counts_per_gene_filtered_hist = px.histogram(gtex_counts_per_gene_filtered, x="x",
									title="GTEx counts per gebe after filtering")
gtex_counts_per_gene_filtered_hist.show()

# TCGA

tcga = pd.read_csv(
	"Data/TCGA/TCGA.tsv", 
	header=0, sep="\t")

tcga_source_site_project = px.histogram(tcga, 
							x="gdc_cases.tissue_source_site.project",
							title="TCGA Number of samples per tissue")
tcga_source_site_project.show()

tcga_age_at_diagnosis = px.histogram(tcga, 
						x="cgc_case_age_at_diagnosis",
						title="TCGA Repartition of ages at diagnosis")
tcga_age_at_diagnosis.show()

tcga_counts = pd.read_csv(
	"Data/TCGA/BeforeNorm/tcga_counts_per_sample.tsv", 
	header=0, sep="\t")

tcga_counts_per_sample = px.histogram(tcga_counts, x="Counts",
							title="TCGA Counts per sample")
tcga_counts_per_sample.show()

tcga_counts_tissue = pd.read_csv(
	"Data/TCGA/BeforeNorm/TCGA_CountsPerTissue.tsv", 
	header=0, sep="\t")

tcga_counts_per_tissue = px.histogram(tcga_counts_tissue, 
							x="Tissue", y='Counts',
							title="TCGA Counts per tissue")
tcga_counts_per_tissue.show()

tcga_mean_counts = pd.read_csv(
	"Data/TCGA/BeforeNorm/tcga_mean_counts_per_sample.tsv", 
	header=0, sep="\t")

tcga_mean_counts_per_sample_hist = px.histogram(tcga_mean_counts, x="Counts",
									title="TCGA Mean counts per sample")
tcga_mean_counts_per_sample_hist.show()

tcga_mean_counts_per_sample_box = px.box(tcga_mean_counts, x="Counts",
									title="TCGA Mean counts per sample")
tcga_mean_counts_per_sample_box.show()

tcga_counts_per_gene = pd.read_csv(
	"Data/TCGA/BeforeNorm/CountsPerGene.tsv",
	header=0, sep="\t")

tcga_counts_per_gene_hist = px.histogram(tcga_counts_per_gene, x="Counts",
							title="TCGA Counts per gene")
tcga_counts_per_gene_hist.show()

tcga_counts_per_gene_box = px.box(tcga_counts_per_gene, x="Counts",
							title="TCGA Counts per gene")
tcga_counts_per_gene_box.show()