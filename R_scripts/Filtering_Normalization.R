library(vsn)
library(edgeR)
library(recount)
library(BatchQC)
library(org.Hs.eg.db)

gtex_raw_counts = load(
	"/Users/labo/Documents/Data/gtex/rdata/PairedEndRounded.Rdata") 

y = DGEList(c1)

# Mean-variance plot
meanSdPlot(y$counts)

# Filter to remove low counts
keep = rowSums(cpm(y) > 5) >= 18  
table(keep)
x = y[keep, ,keep.lib.sizes=FALSE]
write.table(x, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/GTExFilteredCPM5S18.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(x$counts)

# Normalization
x = calcNormFactors(x, method="TMM")
tmm = cpm(x)
write.table(tmm, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/GTExNormalized.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(x$counts)
