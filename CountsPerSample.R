library(recount)

load("/Users/labo/Documents/Data/gtex/rdata/PairedEndRounded.Rdata")

counts_per_samples = colSums(c1)
write.table(counts_per_samples, "/Users/labo/Documents/Data/gtex_counts_per_sample.tsv", sep="\t")

load("/Users/labo/Documents/Data/tcga/rdata/PairedEndRounded.Rdata")

counts_per_samples = colSums(c1)
write.table(counts_per_samples, "/Users/labo/Documents/Data/tcga_counts_per_sample.tsv", sep="\t")