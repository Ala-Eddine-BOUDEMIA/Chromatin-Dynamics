import pandas as pd 
import plotly.express as px

gtex = pd.read_csv(
	"Data/GTEx.tsv", 
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

fig2 = px.histogram(tcga, x="cgc_case_age_at_diagnosis")
fig2.show()