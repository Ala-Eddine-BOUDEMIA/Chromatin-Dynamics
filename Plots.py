import pandas as pd 
import plotly.express as px

gtex = pd.read_csv(
	"Data/GTEx/GTEx.tsv", 
	header=0, sep="\t")

gtex_counts = pd.read_csv(
	"Data/GTEx/gtex_counts_per_sample.tsv", 
	header=0, sep="\t")

gtex_counts_tissue = pd.read_csv(
	"Data/GTEx/GTEx_CountsPerTissue.tsv", 
	header=0, sep="\t")

gtex_mean_counts = pd.read_csv(
	"Data/GTEx/gtex_mean_counts_per_sample.tsv", 
	header=0, sep="\t")

tcga = pd.read_csv(
	"Data/TCGA/TCGA.tsv", 
	header=0, sep="\t")

tcga_counts = pd.read_csv(
	"Data/TCGA/tcga_counts_per_sample.tsv", 
	header=0, sep="\t")

tcga_counts_tissue = pd.read_csv(
	"Data/TCGA/TCGA_CountsPerTissue.tsv", 
	header=0, sep="\t")

tcga_mean_counts = pd.read_csv(
	"Data/TCGA/tcga_mean_counts_per_sample.tsv", 
	header=0, sep="\t")

fig = px.histogram(gtex, x="smtsd")
fig.show()

fig1 = px.histogram(gtex, x="smts")
fig1.show()

fig2 = px.histogram(tcga, x="gdc_cases.tissue_source_site.project")
fig2.show()

fig3 = px.histogram(tcga, x="cgc_case_age_at_diagnosis")
fig3.show()

fig4 = px.histogram(gtex_counts, x="x")
fig4.show()

fig5 = px.box(gtex_counts, x="x")
fig5.show()

fig6 = px.histogram(gtex_mean_counts, x="x")
fig6.show()

fig7 = px.box(gtex_mean_counts, x="x")
fig7.show()

fig8 = px.histogram(gtex_counts_tissue, x="Tissue", y='Counts')
fig8.show()

fig9 = px.histogram(tcga_counts, x="x")
fig9.show()

fig10 = px.box(tcga_counts, x="x")
fig10.show()

fig11 = px.histogram(tcga_counts_tissue, x="Tissue", y='Counts')
fig11.show()

fig12 = px.histogram(tcga_mean_counts, x="x")
fig12.show()

fig12 = px.box(tcga_mean_counts, x="x")
fig12.show()