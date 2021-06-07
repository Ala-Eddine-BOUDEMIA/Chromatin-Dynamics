import pandas as pd 
from itertools import islice

gtex = pd.read_csv(
	"Data/GTEx/GTEx.tsv", 
	header=0, sep="\t")

tcga = pd.read_csv(
	"Data/TCGA/TCGA.tsv", 
	header=0, sep="\t")

gtex_runs_per_tissue = {}
for i, j in zip(gtex["smts"], gtex["run"]):
	if i in gtex_runs_per_tissue.keys():
		gtex_runs_per_tissue[i].append(j)
	else:
		gtex_runs_per_tissue[i] = [j]

with open("Data/GTEx/GTEx_RunsPerTissue.tsv","w") as w:
	for i, j in zip(gtex_runs_per_tissue.keys(), gtex_runs_per_tissue.values()):
		w.write(str(i) + "\t" + str(j) + "\n") 

gtex_counts_per_tissue = {}
with open("Data/GTEx/gtex_counts_per_sample.tsv","r") as f:
	for line in islice(f,1,None):
		l = line.split("\t")
		for i, k in zip(gtex_runs_per_tissue.values(), gtex_runs_per_tissue.keys()):
			for j in i:
				if str(l[0]) == str('"') + j + str('"'):
					if k in gtex_counts_per_tissue.keys(): 
						gtex_counts_per_tissue[k] += int(l[1])
					else:
						gtex_counts_per_tissue[k] = int(l[1])

with open("Data/GTEx/GTEx_CountsPerTissue.tsv","w") as w:
	w.write("Tissue"+ "\t" + "Counts" + "\n")
	for i, j in zip(gtex_counts_per_tissue.keys(), gtex_counts_per_tissue.values()):
		w.write(str(i) + "\t" + str(j) + "\n") 

tcga_runs_per_tissue = {}
for i, j in zip(tcga["gdc_cases.project.primary_site"], tcga["bigwig_file"]):
	if i in tcga_runs_per_tissue.keys():
		tcga_runs_per_tissue[i].append(j)
	else:
		tcga_runs_per_tissue[i] = [j]

with open("Data/TCGA/TCGA_FilesPerTissue.tsv","w") as w:
	for i, j in zip(tcga_runs_per_tissue.keys(), tcga_runs_per_tissue.values()):
		w.write(str(i) + "\t" + str(j) + "\n") 

tcga_counts_per_tissue = {}
with open("Data/TCGA/tcga_counts_per_sample.tsv","r") as f:
	for line in islice(f,1,None):
		l = line.split("\t")
		for i, k in zip(tcga_runs_per_tissue.values(), tcga_runs_per_tissue.keys()):
			for j in i:
				try:
					if str(l[0]) == str('"') + j.split(".")[0] + str('"'):
						if k in tcga_counts_per_tissue.keys(): 
							tcga_counts_per_tissue[k] += int(l[1])
						else:
							tcga_counts_per_tissue[k] = int(l[1])
				except:
					pass

with open("Data/TCGA/TCGA_CountsPerTissue.tsv","w") as w:
	w.write("Tissue"+ "\t" + "Counts" + "\n")
	for i, j in zip(tcga_counts_per_tissue.keys(), tcga_counts_per_tissue.values()):
		w.write(str(i) + "\t" + str(j) + "\n") 