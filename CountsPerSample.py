import pandas as pd 

gtex_counts = pd.read_csv(
	"Data/GTEx_counts.tsv", 
	header=0, index_col=0, sep="\t")

d = {}
for i in gtex_counts.columns:
	s = 0
	d[i] = 0
	for j in gtex_counts[i]:
		s += int(j)
	d[i] = s

with open("Data/countsPerSample.tsv","w") as w:
	for i,j in zip(d.keys(), d.values()):
		w.write(str(i) + "\t" + str(j) + "\n") 