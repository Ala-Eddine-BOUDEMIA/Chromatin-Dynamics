library(recount)
tcga = load("/Users/labo/Documents/Data/rse_gene_TCGA.Rdata")
gtex = load("/Users/labo/Documents/Data/rse_gene_GTEx.Rdata")
tcga_counts = assays(read_counts(tcga, TRUE, TRUE))$counts
gtex_counts = assays(read_counts(gtex, TRUE, TRUE))$counts