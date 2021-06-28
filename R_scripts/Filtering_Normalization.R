library(vsn)
library(edgeR)
library(recount)
library(BatchQC)
library(org.Hs.eg.db)

gtex_raw_counts = load(
	"/Users/labo/Documents/Data/gtex/rdata/PairedEndRounded.Rdata") 

gtex_metadata = read.delim(
	"/Users/labo/Documents/Data/gtex/tsv/GTEx.tsv", 
	header=TRUE, sep="\t")

tissues = gtex_metadata$smts
tissues = factor(tissues)

y = DGEList(c1, group=tissues, genes=c1[,1,drop=FALSE])

y$genes$Symbol = mapIds(org.Hs.eg.db,
                        keys = rownames(y), 
                        keytype = 'ENSEMBL', 
                        column = 'SYMBOL')
head(y$genes)
y = y[!is.na(y$genes$Symbol), ]
dim(y)

# Mean-variance plot
meanSdPlot(y$counts)

# Filter to remove low counts
keep = rowSums(cpm(y) > 10) >= 36  
table(keep)
x = y[keep, ,keep.lib.sizes=FALSE]
write.table(x, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/GTExFilteredCPM10S36.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(x$counts)

# Normalization
x = calcNormFactors(x)
x$samples
write.table(x$counts, "/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/GTExNormalized.tsv", sep="\t")

# Mean-variance plot
meanSdPlot(x$counts)
