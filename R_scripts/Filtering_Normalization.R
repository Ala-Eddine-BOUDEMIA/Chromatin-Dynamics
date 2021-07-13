library(vsn)
library(edgeR)
library(recount)
library(sva)
library(org.Hs.eg.db)

gtex_filtered = read.delim("/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/ABER/Counts/AfterFiltering/Filtered_smnabtch.tsv", 
                           header = TRUE, sep = '\t')

rownames(gtex_filtered) = gtex_filtered$X

metadata = read.delim('/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/ABER/Metadata/GTEx_smnabtch.tsv',
                      header = TRUE, sep = "\t")

rownames(metadata) = metadata$run

y = DGEList(data.matrix(gtex_filtered))

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
batch = metadata$smnabtch
modcombat = model.matrix(~1, data=batch)
gtex_filtered = as.matrix(gtex_filtered)
corrected_counts = ComBat_seq(gtex_filtered, batch=batch)

v = DGEList(data.matrix(corrected_counts))

# Normalization
z = calcNormFactors(v, method="TMM")
tmm = cpm(z)
write.table(tmm, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/ABER/Counts/Normalized/NormalizedCorrected.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(tmm)
