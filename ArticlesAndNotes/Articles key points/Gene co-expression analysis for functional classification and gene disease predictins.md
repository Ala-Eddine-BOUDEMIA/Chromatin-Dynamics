# Gene co-expression analysis for functional classification and gene–disease predictions

Citation: van Dam, S., Võsa, U., van der Graaf, A., Franke, L., & de Magalhães, J. P. (2018). Gene co-expression analysis for functional classification and gene-disease predictions. Briefings in bioinformatics, 19(4), 575–592. https://doi.org/10.1093/bib/bbw139

- Gene co-expression networks can be used for various purposes, including candidate disease gene prioritization, functional gene annotation, and the identification of regulatory genes. However, co-expression networks are effectively only able to identify correlations; they indicate which genes are active simultaneously, which often indicates they are active in the same biological processes but do not normally confer information about causality or distinguish between regulatory and regulated genes.
- An increasingly used method that goes beyond traditional co-expression networks is differential co-expression analysis. This approach identifies genes with varying co-expression partners under different conditions, such as disease states, tissue types, and developmental stages, because these genes are more likely to be regulators that underlie phenotypic differences.
- Gene expression and regulation can be highly tissue-specific, and most disease-related genes have tissue-specific expression abnormalities. The increased availability of expression data for multiple tissues has allowed for differential co-expression analysis, which can identify both tissue-specific signatures and shared co-expression signatures.

# Co-expression networks

- A co-expression network identifies which genes have a tendency to show a coordinated expression pattern across a group of samples.
- In the first step, individual relationships between genes are defined based on correlation measures or mutual information between each pair of genes. These relationships describe the similarity between expression patterns of the gene pair across all the samples.
- In the second step, co-expression associations are used to construct a network where each node represents a gene and each edge represents the presence and the strength of the co-expression relationship
- In the third step, modules (groups of co-expressed genes) are identified using one of several available clustering techniques.
- In co-expression analysis, it is important to consider the heterogeneity of the samples. Tissue-specific or condition-specific co-expression modules may not be detectable in a co-expression network constructed from multiple tissues or conditions because the correlation signal of the tissue/ condition-specific modules is diluted by a lack of correlation in other tissues/conditions. However, limiting co-expression analysis to a specific tissue or condition also reduces the sample size, thereby also decreasing the statistical power to detect shared co-expression modules. Therefore, methods that do not distinguish between tissues or conditions should be used for the identification of common co-expression modules, while differential co-expression comparing different conditions or tissues will be better for identifying modules unique to a specific condition or tissue.

## Types of co-expression networks

### Signed and unsigned co-expression networks

### Weighted and unweighted co-expression networks

## Microarrays Vs. RNA-seq data

- One of the major benefits of RNA-seq is that it quantifies the expression of the over 70 000 non-coding RNAs not usually measured with microarrays
- RNA-seq also has other benefits. It increases accuracy for low-abundance transcripts, has a higher resolution for identifying tissue-specific expression, and distinguishes expression profiles of closely related paralogues better than microarray-derived profiles. RNA-seq can also distinguish between the expression of different splice variants, which can have distinct interaction partners and biological functions.

# Clustering and network analysis

## Identifying modules

- Clustering is used to group genes that have a similar expression pattern in multiple samples. The resulting modules often represent biological processes and can be phenotype-specific.
- Hierarchical clustering iteratively divides each cluster into sub-clusters to create a tree with branches representing co-expression modules. Modules are then defined by cutting the branches at a certain height

## Identifying hub genes

- Hub genes are highly connected genes in a co-expression network.
- Hubs are frequently more relevant to the functionality of networks than other nodes.

## Guilt by association

- A widely used approach to attach biological meaning to modules is to determine functional enrichment among the genes within a module.
- Assuming that co-expressed genes are functionally related, enriched functions can be assigned to poorly annotated genes within the same co-expression module, an approach commonly referred to as ‘guilt by association.

# Differential co-expression analysis

- Genes that are differentially co-expressed between different sample groups are more likely to be regulators and are therefore likely to explain differences between phenotypes

## Differential co-expression analysis between sample groups

- Most differential co-expression analyses rely on differential clustering; they identify clusters that contain different genes or behave differently under changing conditions or phenotypes.
- Programs for differential co-expression identify co-expressed modules across the full set of samples.
- The co-expressed modules can then be correlated to predefined sample sub-populations (disease or tissue ...etc).
- WGCNA determines the activity and importance of each module in each sub-population of samples :
    - For each module, an eigengene is calculated, which is the vector that best describes the expression behavior of all genes within this module.
    - It then prioritizes which genes in these modules are likely to underlie the phenotype associated with the module by :
        - Identifying which genes behave similarly to the eigengene of the module
        - Identifying which genes that are intra-modular hub genes
- DICER is tailored to identify module pairs that correlate differently between sample groups.
- DICER is useful for time series experiments in which co-expression changes are gradual.
- DiffCoEx focuses on modules that are differentially co-expressed with the same sets of genes.