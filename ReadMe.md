# M1 GENIOMHE Internship at Chromatin-Dynamics Team
![Solid visualization](/ReadMe_Images/Logos.png)

**Integrated Analysis of Co-expression Patterns of Histone Variants and their Chaperones Across Tissues from Publicly Available RNA-Seq Data**


# Table of contents

* [General Information](#General-Information)
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Usage](#Usage)

# General Information

## Aim of the project

Biochemical studies have allowed the characterization of several interactions between histone variants and their chaperones. However, their regulation across and within different tissues is still not fully characterized neither in healthy nor in disease context.
The objective of this project was to perform a co-expression analysis using publicly accessible transcriptomics data to investigate how histone chaperones and histone variants are expressed across and within tissues and how these expression patterns would change in the context of distinct cancers.

## Data Acquisition

The data was downloaded from [recount2 website](https://jhubiostatistics.shinyapps.io/recount/).
The Two datasets used in this analysis are the Genotype-Tissue Expression (GTEx) Dataset and The Cancer Genome Atlas (TCGA) Dataset.
The data provided by recount2 is per-base coverage, it was transformed into normal counts using the script [init_data.R](/R_scripts/init_data.R).
The data should be organized as follow:

    Chromatin-Dynamics
                    └─── Data
                            └─── dataset (GTEx or TCGA) 				
                                                    └──── BeforeFiltering
                                                                        └──── PairedEndRounded.tsv


# Requirements

- dash==1.20.0
- dash_bio==0.7.0
- matplotlib==3.4.2
- networkx==2.5.1
- numpy==1.21.0
- pandas==1.2.4
- plotly==4.14.3
- pyvis==0.1.9
- scikit-learn==0.24.2
- scipy==1.6.3
- seaborn==0.11.1
- umap==0.1.1

# Installation 

## Install packages:


```
sudo pip3 install -r requirements.txt
```


## Download repository:

```
git clone https://github.com/Ala-Eddine-BOUDEMIA/Chromatin-Dynamics.git
```


```
cd Chromatin-Dynamics
```

# Usage

- Take a look at `Config.py` before you begin to get a feel for what parameters can be changed.
- Make sure to update the parameter `dataset` depending on the dataset you want to use

## CPM

- The first thing to do is to normalize the data to the sequencing depth, therfore run the following command.

```
python3 CPM.py
```

**Inputs**: `Data/dataset/BeforeFiltering/PairedEndRounded.tsv` where dataset == GTEx or TCGA

**Outputs**: `Data/dataset/CPM/Counts/Normalized/Full/counts.tsv`

- Note that there is a filtering command to filter lowly expressed genes, feel free to change the threshold.
- For TCGA it will also filter the samples; possible sample types are:
  - Primary Tumor
  - Recurrent Tumor
  - Metastatic
  - Solid Tissue Normal
  - Primary Blood Derived Cancer - Peripheral Blood

    Chromatin-Dynamics
        └─── Data
                └─── dataset (GTEx or TCGA) 			
                        └──── CPM				
                                └──── Counts			
                                        └──── Normalized		
                                                └──── Full   
                                                        └──── Counts.tsv

## 0_GenerateData.py

- This code uses the normalized counts to generate subsections of the datasets.

**Inputs**: `Data/GTEx/CPM/Counts/Normalized/Full/counts.tsv` 

If the parameter `which` in the `Config.py` file is set to **Normalized** 

**Outputs**: `Data/GTEx/CPM/Counts/Normalized/CountsByTissue/**/*.tsv`, `Data/GTEx/CPM/Counts/Normalized/Normal/counts.tsv`, `Data/GTEx/CPM/Counts/Normalized/WoTissues/counts.tsv`, `Data/GTEx/CPM/Counts/Top1000/counts.tsv`, `Data/GTEx/CPM/Counts/Top100/counts.tsv`, `Data/GTEx/CPM/Counts/Random/*.tsv`

- CountsByTissue: Contains separate files for each tissue type with the counts for all the genes 
- Normal: Drops the samples coming from the transformed cells
- WoTissues: Drops the samples coming from one of these tissues: Brain, Blood, Bone Marrow, Pituitary, Spleen, Testis
- Top1000: Contains all the samples but only the top 1000 expressed genes
- Top100: Contains all the samples but only the top 125 expressed genes
- Random: Contains 10 files each contain all the samples and only 125 randomly selected genes

If the parameter `which` in the `Config.py` file is set to **variants_chaperones** the outputs for GTEx will be: 

**Outputs**: `Data/GTEx/CPM/Counts/variants_chaperones/Full/counts.tsv`, `Data/GTEx/CPM/Counts/variants_chaperones/NonReplicative/counts.tsv`, `Data/GTEx/CPM/Counts/variants_chaperones/CountsByTissue/**/*.tsv`, `Data/GTEx/CPM/Counts/variants_chaperones/Normal/counts.tsv`, `Data/GTEx/CPM/Counts/variants_chaperones/WoTissues/counts.tsv`, `Data/GTEx/CPM/Counts/Top1000/counts.tsv`, `Data/GTEx/CPM/Counts/Top100/counts.tsv`, `Data/GTEx/CPM/Counts/Random/*.tsv`

- Full: Contains all the counts for all the histone chaperone and histone variant genes
- NonReplicative: Contains all the counts for all the histone chaperone and non-replicative histone variant genes
- CountsByTissue: Contains separate files for each tissue type with counts of the histone chaperone and histone variant genes only
- Normal: Drops the samples coming from the transformed cells
- WoTissues: Drops the samples coming from one of these tissues: Brain, Blood, Bone Marrow, Pituitary, Spleen, Testis
- Top1000: Contains all the samples but only the top 1000 expressed genes
- Top100: Contains all the samples but only the top 125 expressed genes
- Random: Contains 10 files each contain all the samples and only 125 randomly selected genes

**Inputs**: `Data/TCGA/CPM/Counts/Normalized/Full/counts.tsv` 

If the parameter `which` in the `Config.py` file is set to **Normalized** 

**Outputs**: `Data/TCGA/CPM/Counts/Normalized/CountsByTissue/**/*.tsv`, `Data/TCGA/CPM/Counts/Top1000/counts.tsv`, `Data/TCGA/CPM/Counts/Top100/counts.tsv`, `Data/TCGA/CPM/Counts/Random/*.tsv`

- CountsByTissue: Contains separate files for each cancer type with the counts for all the genes 
- Top1000: Contains all the samples but only the top 1000 expressed genes
- Top100: Contains all the samples but only the top 125 expressed genes
- Random: Contains 10 files each contain all the samples and only 125 randomly selected genes

If the parameter `which` in the `Config.py` file is set to **variants_chaperones** the outputs for GTEx will be: 

**Outputs**: `Data/TCGA/CPM/Counts/variants_chaperones/Full/counts.tsv`, `Data/TCGA/CPM/Counts/variants_chaperones/NonReplicative/counts.tsv`, `Data/TCGA/CPM/Counts/variants_chaperones/CountsByTissue/**/*.tsv`, `Data/TCGA/CPM/Counts/Top1000/counts.tsv`, `Data/TCGA/CPM/Counts/Top100/counts.tsv`, `Data/TCGA/CPM/Counts/Random/*.tsv`

- Full: Contains all the counts for all the histone chaperone and histone variant genes
- NonReplicative: Contains all the counts for all the histone chaperone and non-replicative histone variant genes
- CountsByTissue: Contains separate files for each cancer type with counts of the histone chaperone and histone variant genes only
- Top1000: Contains all the samples but only the top 1000 expressed genes
- Top100: Contains all the samples but only the top 125 expressed genes
- Random: Contains 10 files each contain all the samples and only 125 randomly selected genes

   Chromatin-Dynamics
        └─── Data
                └─── dataset (GTEx or TCGA) 			
                        └──── CPM				
                                └──── Counts			
                                        └──── Normalized or Random or Top1000 or Top100 or variants_chaperones	
                                                └──── Normal or WoTissues or CountsByTissue    
                                                        └──── Counts.tsv

```
python3 0_GenerateData.py 
```

## 1_ExploreData.py

- This code generates plots to see:
  1. How the counts per sample are distributed
  2. How the counts per gene are distributed
  3. How the counts per sample are distributed within each tissue/cancer type
  4. How the counts per gene are distributed within each tissue/cancer type

- As for the previous code make sure to update the parameter `which` in the `Config.py` file depending on what set of genes you want to work on.
- Also make sure to select which subsection of the dataset you want to execute the code on by updating the `counts`, `general_qc_imgs` and `general_qc_htmls` in the `Config.py` file

**Inputs**: `counts`
**Outputs**: `general_qc_imgs` and `general_qc_htmls`

```
python3 1_ExploreData.py 
```

## 2_MeanStd.py

- This code generate the mean_variance plots on the log2(CPMs + 1)
- This code is mainly useful on the complete set of genes therefore it is recomended to run it only with the parameter `which` set to **Normalized**
- Make sure to select which subsection of the dataset you want to execute the code on by updating the `counts`, `mv_imgs` and `mv_htmls` in the `Config.py` file

**Inputs**: `counts`
**Outputs**: `mv_imgs` and `mv_htmls`

```
python3 2_MeanStd.py 
```

## 3_PCA.py

- This code generates PCA plots and colors the samples by their tissue type.
- Make sure to update the parameter `which` in the `Config.py` file depending on what set of genes you want to work on.
- Make sure to select which subsection of the dataset you want to execute the code on by updating the `counts`, `pca_imgs`, `pca_htmls` and `files_pca` in the `Config.py` file

**Inputs**: `counts`
**Outputs**: `pca_imgs`, `pca_htmls` and `files_pca` 

```
python3 3_PCA.py 
```

## 4_TSNE.py

- This code generates t-SNE plots and colors the samples by their tissue type.
- Make sure to update the parameter `which` in the `Config.py` file depending on what set of genes you want to work on.
- Make sure to select which subsection of the dataset you want to execute the code on by updating the `counts`, `tsne_imgs`, `tsne_htmls` and `files_tsne` in the `Config.py` file

**Inputs**: `counts`
**Outputs**: `tsne_imgs`, `tsne_htmls` and `files_tsne` 

```
python3 4_TSNE.py 
```

## 5_UMAP.py

- This code generates UMAP plots and colors the samples by their tissue type.
- Make sure to update the parameter `which` in the `Config.py` file depending on what set of genes you want to work on.
- Make sure to select which subsection of the dataset you want to execute the code on by updating the `counts`, `umap_imgs`, `umap_htmls` and `files_umap` in the `Config.py` file

**Inputs**: `counts`
**Outputs**: `umap_imgs`, `umap_htmls` and `files_umap` 

```
python3 5_UMAP.py 
```

## 6_Correlation.py

- This code generates Pearson correlation matrices.
- Make sure to update the parameter `which` in the `Config.py` file depending on what set of genes you want to work on.
- Make sure to select which subsection of the dataset you want to execute the code on by updating the `counts`, `s_corr` and `g_corr` in the `Config.py` file

**Inputs**: `counts`
**Outputs**: `s_corr` and `g_corr`

- For files that does not contain enough samples it is possible that generated files will contain *NaN* values make sure to eliminate them.

```
python3 6_Correlation.py 
```

## 7_Clustering.py

- This code generates sample by sample, gene by gene and sample by gene clustermaps.
- Make sure to update the parameter `which` in the `Config.py` file to **variants_chaperones** because the **Normalized** set is too huge for this computation.
- Make sure to select which subsection of the dataset you want to execute the code on by updating the `counts`, `s_corr`, `g_corr`, `s_clustermaps`, `g_clustermaps` and `s_clustermaps` in the `Config.py` file

**Inputs**: `counts`, `s_corr` and `g_corr`
**Outputs**: `s_clustermaps`, `g_clustermaps` and `s_clustermaps`

```
python3 7_Clustering.py 
```

- Note that there is a variable called `c`, it should be modified depending on the dataset you enter otherwise you might get an error or a clustermap without gene names.

## 8_AdjacencyMatrix.py

- This code takes an adjacency matrix and generates an interactive graph
- The generated graphs could be found at `WGCNA/dataset/Networks/*.html`
