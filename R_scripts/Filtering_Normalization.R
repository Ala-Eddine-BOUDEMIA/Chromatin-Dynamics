library(vsn)
library(edgeR)
library(recount)
library(sva)
library(org.Hs.eg.db)

gtex_raw_counts = load(
	"/Users/labo/Documents/Data/gtex/rdata/PairedEndRounded.Rdata")

metadata = read.delim('/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/BBER/Metadata/GTEx.tsv',
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

# Batch Effect Removal
batch1 = metadata$smnabtch
batch2 = metadata$smgebtch

corrected_batch1 = ComBat(x, batch1)
corrected_batch1_2 = ComBat(corrected_batch1, batch2)

# Normalization
z = calcNormFactors(combat_data2, method="TMM")
tmm = cpm(z)
write.table(tmm, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/GTExNormalizedCorrected.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(tmm)
