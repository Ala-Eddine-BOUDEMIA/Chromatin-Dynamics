library(vsn)
library(edgeR)
library(recount)
library(sva)
library(org.Hs.eg.db)

gtex_raw_counts = load(
	"/Users/labo/Documents/Data/gtex/rdata/PairedEndRounded.Rdata")

metadata = read.delim('/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/Metadata/GTEx.tsv',
                      header = TRUE, sep = "\t")

rownames(metadata) = metadata$run

y = DGEList(c1)

# Mean-variance plot
meanSdPlot(y$counts)

# Filter to remove low counts
keep = rowSums(cpm(y) > 5) >= 18  
table(keep)
x = y[keep, ,keep.lib.sizes=FALSE]
write.table(x$counts, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/GTExFilteredCPM5S18.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(x$counts)

# Normalization
z = calcNormFactors(combat_data2, method="TMM")
tmm = cpm(z)
write.table(tmm, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/GTExNormalizedCorrected.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(tmm)
