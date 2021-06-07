import pandas as pd 
import plotly.express as px

gtex = pd.read_csv(
	"Data/GTEx.tsv", 
	header=0, sep="\t")

gtex_counts = pd.read_csv(
	"Data/gtex_counts_per_sample.tsv", 
	header=0, sep="\t")

gtex_counts_tissue = pd.read_csv(
	"Data/CountsPerTissue.tsv", 
	header=0, sep="\t")

tcga = pd.read_csv(
	"Data/TCGA.tsv", 
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

fig6 = px.histogram(gtex_counts_tissue, x="Tissue", y='Counts')
fig6.show()