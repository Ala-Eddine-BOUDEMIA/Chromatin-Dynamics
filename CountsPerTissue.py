import pandas as pd 
from itertools import islice

gtex_counts = pd.read_csv(
	"Data/GTEx.tsv", 
	header=0, sep="\t")

runs_per_tissue = {}
for i, j in zip(gtex_counts["smts"], gtex_counts["run"]):
	if i in runs_per_tissue.keys():
		runs_per_tissue[i].append(j)
	else:
		runs_per_tissue[i] = [j]

with open("Data/RunsPerTissue.tsv","w") as w:
	for i, j in zip(runs_per_tissue.keys(), runs_per_tissue.values()):
		w.write(str(i) + "\t" + str(j) + "\n") 

counts_per_tissue = {}
with open("Data/gtex_counts_per_sample.tsv","r") as f:
	for line in islice(f,1,None):
		l = line.split("\t")
		for i, k in zip(runs_per_tissue.values(), runs_per_tissue.keys()):
			for j in i:
				if str(l[0]) == str('"') + j + str('"'):
					if k in counts_per_tissue.keys(): 
						counts_per_tissue[k] += int(l[1])
					else:
						counts_per_tissue[k] = int(l[1])

with open("Data/CountsPerTissue.tsv","w") as w:
	w.write("Tissue"+ "\t" + "Counts" + "\n")
	for i, j in zip(counts_per_tissue.keys(), counts_per_tissue.values()):
		w.write(str(i) + "\t" + str(j) + "\n") 