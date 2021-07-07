library(BatchQC)
library(sva)

metadata = read.delim('/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/Metadata/GTEx.tsv',
                      header = TRUE, sep = "\t")

rownames(metadata) = metadata$run

normalized_counts = read.delim('/Users/labo/Documents/Code/Chromatin-Dynamics/Data/GTEx/Counts/Normalized/Normalized.tsv',
                               header = TRUE, sep = "\t") 

nsample = 9662
sample = colnames(normalized_counts)
condition = metadata$smts
batch = data.frame(metadata$smnabtch, metadata$smgebtch, metadata$smtspax, metadata$smtstptref)
rownames(batch) = metadata$run

pdata <- data.frame(sample, batch, condition)
modmatrix = model.matrix(~as.factor(condition), data=pdata)
combat_data.matrix = ComBat(normalized_counts, batch=metadata$smnabtch)
