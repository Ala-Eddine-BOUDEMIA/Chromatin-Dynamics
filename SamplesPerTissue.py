import pandas as pd 
import plotly.express as px

gtex = pd.read_csv(
	"Data/GTEx.tsv", 
	header=0, sep="\t")

fig = px.histogram(gtex["smtsd"], x="smtsd")
fig.show()