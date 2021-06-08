library(edgeR)
library("org.Hs.eg.db")

gtex_raw_counts = read.delim(
	"/Users/labo/Documents/Data/gtex/tsv/PairedEndRounded.tsv", 
	header=TRUE, sep="\t") 
	#row.names="EnsembleGeneID" see how it is written

gtex_metadata = read.delim(
	"/Data/GTEx/GTEx.tsv", 
	header=TRUE, sep="\t")

tissues = gtex_metadata$smtsd
tissues = factor(tissues)

y = DGEList(gtex_raw_counts[,-1],
		group=tissues,
		genes=gtex_raw_counts[,1,drop=FALSE]) #check the output of this line

y$genes$Symbol = mapIds(org.Mm.eg.db, 
					rownames(y), 
					keytype="ENSEMBL", 
					column="SYMBOL")
head(y$genes)
y = y[!is.na(y$genes$Symbol), ]
dim(y)

# Filter to remove low counts
# try different cpms
# https://f1000researchdata.s3.amazonaws.com/manuscripts/9996/f522df8c-8041-48a8-afb1-038531488f95_8987_-_gordon_smyth_v2.pdf?doi=10.12688/f1000research.8987.2&numberOfBrowsableCollections=28&numberOfBrowsableInstitutionalCollections=4&numberOfBrowsableGateways=25
# CPM has been chosen because it is roughly equal to 10/L where L is the minimum library size in millions. 
# In our case 10/3.4 = 3
keep = rowSums(cpm(y) > 3) >= 2 
table(keep)
y[keep, ,keep.lib.sizes=FALSE]

# QC1
## multi-dimensional scaling
#pch = c(0,1,2,15,16,17)
#colors = rep(c("darkgreen","red","blue"),2)
plotMDS(y)#,col=colors[group],pch=pch[group])
#legend("topleft",legend=levels(group),pch=pch,col=colors,ncol=2)

## mean-difference
#It is good practice to make MD plots for all the samples as a quality check. column 1 to +9000
plotMD(y,column=1)
abline(h=0,col="red",lty=2,lwd=2)

# Normalization
y = calcNormFactors(y)
y$samples

# QC2
## multi-dimensional scaling
#pch = c(0,1,2,15,16,17)
#colors = rep(c("darkgreen","red","blue"),2)
plotMDS(y)#,col=colors[group],pch=pch[group])
#legend("topleft",legend=levels(group),pch=pch,col=colors,ncol=2)

## mean-difference
#It is good practice to make MD plots for all the samples as a quality check. column 1 to +9000
plotMD(y,column=1)
abline(h=0,col="red",lty=2,lwd=2)