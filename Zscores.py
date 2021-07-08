import pandas as pd
import plotly.express as px

import Config

from sklearn.preprocessing import StandardScaler

counts = pd.read_csv(Config.args.cv, 
	header = 0, index_col = 0, sep = "\t")

cv = pd.read_csv(Config.args.list, 
	header = 0, index_col = 0, sep = ";")

scaler = StandardScaler()
std_counts = pd.DataFrame(scaler.fit_transform(counts), 
	index = counts.index, columns = counts.columns)

df = pd.DataFrame(index = counts.index)
mean = pd.DataFrame(std_counts.mean(axis = 1), columns = ["mean"])
std = pd.DataFrame(std_counts.std(axis = 1), columns = ['std']) 

df = df.join(mean)
df = df.join(std)

df = df.sort_values("mean")
#df = df.join(cv["GeneName"])

fig = px.histogram(df)
fig.show()