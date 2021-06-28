library(recount)
tcga = load("/Users/labo/Documents/tcga/rdata/Data/rse_gene_TCGA.Rdata")
gtex = load("/Users/labo/Documents/gtex/rdata/Data/rse_gene_GTEx.Rdata")
tcga_counts = assays(read_counts(tcga, TRUE, TRUE))$counts
gtex_counts = assays(read_counts(gtex, TRUE, TRUE))$counts
write.table(tcga_counts, file="/Users/labo/Documents/Data/tcga/PairedEndRounded.tsv", sep="\t")
write.table(gtex_counts, file="/Users/labo/Documents/Data/gtex/PairedEndRounded.tsv", sep="\t")