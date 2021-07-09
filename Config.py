import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

# GTEx or TCGA
parser.add_argument("--dataset",
	type = str,
	default = "GTEx",
	help = "Dataset to use GTEx or TCGA")

# Before or After Batch Effect Removal (A/BBER)
parser.add_argument("--batch",
	type = str,
	default = "BBER",
	help = "Dataset to use before or after batch effect removal")

# Data 
choice = parser.parse_args()

## Counts
parser.add_argument("--bf",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/BeforeFiltering/PairedEndRounded.tsv"),
	help = "Raw counts")

parser.add_argument("--af",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/AfterFiltering/FilteredCPM5S18.tsv"),
	help = "Filtered counts, you can change the file name to one of these: Filtered + CPM5S18 or CPM10S18 or CPM10S36 + .tsv")

parser.add_argument("--norm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/Normalized/Normalized.tsv"),
	help = "Normalized counts after filtering, these are the counts used to generate the following count files")

parser.add_argument("--onlyNormal",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/Normal/counts.tsv"),
	help = "Counts missing the samples coming from transformed cells")

parser.add_argument("--WoTissues",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/WithoutTissues/counts.tsv"),
	help = "Counts missing samples coming from Blood, Brain, Bone, Pituitary and Spleen tissues")

parser.add_argument("--top1000",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/Top1000/Top1000Genes.tsv"),
	help = "Counts limited to the top 1000 expressed genes")

parser.add_argument("--top88",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/Top88/Top88Genes.tsv"),
	help = "Counts limited to the top 88 expressed genes")

parser.add_argument("--cv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/variants_chaperones/variants_chaperones_counts.tsv"),
	help = "Counts limited to the histone chaperones and histone variants genes")

parser.add_argument("--nonRcv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/variants_chaperones/nonRvariants_chaperones_counts.tsv"),
	help = "Counts limited to the histone chaperones and histone non-replicative variants genes")

parser.add_argument("--tissue",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/CountsByTissue/"),
	help = "Each tissue's counts are stored in a different file")

parser.add_argument("--rand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/Counts/Random/"),
	help = "10 randomly generated count files that are limited to the top 88 expressed genes")

## Metadata
parser.add_argument("--meta",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Metadata/" + choice.dataset + ".tsv"),
	help = "Metadata file")

parser.add_argument("--list",
	type = Path,
	default = Path("Data/variants_chaperones/complete_list.csv"),
	help = "List of the names and IDs of histone chaperones and histone variants genes")

parser.add_argument("--nonReplicative",
	type = Path,
	default = Path("Data/variants_chaperones/withoutR.txt"),
	help = "List of the names and IDs of histone chaperones and histone non-replicative variants genes")

## Correlation
parser.add_argument("--corrTopG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Genes/Top1000/corr_matrix.tsv"),
	help = "Correlation matrix between the top 1000 expressed genes")

parser.add_argument("--corr88G",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Genes/Top88/corr_matrix.tsv"),
	help = "Correlation matrix between the top 88 expressed genes")

parser.add_argument("--corrCVg",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Genes/variants_chaperones/corr_matrix.tsv"),
	help = "Correlation matrix between the histone chaperones and histone variants genes")

parser.add_argument("--corrNonRcvG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Genes/nonReplicative/corr_matrix.tsv"),
	help = "Correlation matrix between the histone chaperones and histone non-replicative variants genes")

parser.add_argument("--corrTissueG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Genes/CorrByTissue/"),
	help = "Correlation matrices for each tissue between every gene")

parser.add_argument("--corrRandG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Genes/Random/"),
	help = "Correlation matrices between the randomly selected 88 genes")

parser.add_argument("--corrTopS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Samples/Top1000/corr_matrix.tsv"),
	help = "Correlation matrix between all the samples using counts limited to the top 1000 expressed genes")

parser.add_argument("--corr88S",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Samples/Top88/corr_matrix.tsv"),
	help = "Correlation matrix between all the samples using counts limited to the top 88 expressed genes")

parser.add_argument("--corrCVs",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Samples/variants_chaperones/corr_matrix.tsv"),
	help = "Correlation matrix between all the samples using counts limited to the histone chaperones and histone variants genes")

parser.add_argument("--corrNonRcvS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Samples/nonReplicative/corr_matrix.tsv"),
	help = "Correlation matrix between all the samples using counts limited to the histone chaperones and histone non-replicative variants genes")

parser.add_argument("--corrTissueS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Samples/CorrByTissue/"),
	help = "Correlation matrices between all the samples using counts limited to each tissue")

parser.add_argument("--corrRandS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/CorrelationMatrix/Samples/Random/"),
	help = "Correlation matrices between all the samples using the randomly generated files")

## PCA
parser.add_argument("--pcaRaw",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/BeforeFiltering/pca.tsv"),
	help = "Location where the pca matrix of the raw counts is stored")

parser.add_argument("--pcaFiltered",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/AfterFiltering/pca.tsv"),
	help = "Location where the pca matrix of the filtered non normalized counts is stored")

parser.add_argument("--pcaNorm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/Normalized/pca.tsv"),
	help = "Location where the pca matrix of the filtered normalized counts is stored")

parser.add_argument("--pcaNormal",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/Normal/pca.tsv"),
	help = "Location where the pca matrix of the dataset without transformed cells is stored")

parser.add_argument("--pcaWoTissues",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/WithoutTissues/pca.tsv"),
	help = "Location where the pca matrix for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--pcaTop",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/Top1000/pca.tsv"),
	help = "Location where the pca matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--pcaTop88",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/Top88/pca.tsv"),
	help = "Location where the pca matrix of the top 88 expressed gene counts is stored")

parser.add_argument("--pcaCV",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/variants_chaperones/pca.tsv"),
	help = "Location where the pca matrix of the variant and chaperone gene counts is stored")

parser.add_argument("--pcaTissue",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/PCAbyTissue/"),
	help = "Location where the pca matrix of the gene counts by tissue are stored")

parser.add_argument("--pcaRand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/PCA/Random/"),
	help = "Location where the pca matrix of the randomly selected geen counts is stored")

## T-SNE
parser.add_argument("--tsneRaw",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/BeforeFiltering/T-Sne.tsv"),
	help = "Location where the T-Sne matrix of the raw counts is stored")

parser.add_argument("--tsneFiltered",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/AfterFiltering/T-Sne.tsv"),
	help = "Location where the T-Sne matrix of the filtered non normalized counts is stored")

parser.add_argument("--tsneNorm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/Normalized/T-Sne.tsv"),
	help = "Location where the t-sne matrix of the filtered normalized counts is stored")

parser.add_argument("--tsneNormal",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/Normal/tsne.tsv"),
	help = "Location where the t-sne matrix of the dataset without transformed cells is stored")

parser.add_argument("--tsneWoTissues",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/WithoutTissues/tsne.tsv"),
	help = "Location where the t-sne matrix for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--tsneTop",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/Top1000/T-Sne.tsv"),
	help = "Location where the t-sne matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--tsneTop88",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/Top88/T-Sne.tsv"),
	help = "Location where the t-sne matrix of the top 88 expressed gene counts is stored")

parser.add_argument("--tsneCV",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/variants_chaperones/T-Sne.tsv"),
	help = "Location where the t-sne matrix of the variant and chaperone gene counts is stored")

parser.add_argument("--tsneTissue",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T-Sne/TSNEbyTissue/"),
	help = "Location where the t-sne matrix of the gene counts by tissue are stored")

parser.add_argument("--tsneRand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + choice.batch + "/T_Sne/Random/"),
	help = "Location where the T-Sne matrix of the randomly selected geen counts is stored")

# Images
## General
parser.add_argument("--IgeneralRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/BeforeFiltering/"),
	help = "Location where QC images for the raw dataset are stored")

parser.add_argument("--IgeneralFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/AfterFiltering/"),
	help = "Location where QC images for the filtered dataset are stored")

parser.add_argument("--IgeneralNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/Normalized/"),
	help = "Location where QC images for the normalized dataset are stored")

parser.add_argument("--IgeneralNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/Normal/"),
	help = "Location where QC images for the dataset without transformed cells are stored")

parser.add_argument("--IgeneralWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/WithoutTissues/"),
	help = "Location where QC images for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--IgeneralTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/Top1000/"),
	help = "Location where QC images for the top1000 expressed genes dataset are stored")

parser.add_argument("--IgeneralTop88",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/Top88/"),
	help = "Location where QC images for the top 88 expressed genes dataset are stored")

parser.add_argument("--IgeneralCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/variants_chaperones/"),
	help = "Location where QC images for the chaperone and variants dataset are stored")

parser.add_argument("--IgeneralNrCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/NonReplicativeHistones/"),
	help = "Location where QC images for the chaperone and non replicative variants dataset are stored")

parser.add_argument("--IgeneralRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/General/Random/"),
	help = "Location where QC images for the Random dataset are stored")

## Mean Variance
parser.add_argument("--ImvRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/BeforeFiltering/mv.png"),
	help = "Location where the mean-variance image of the raw dataset is stored")

parser.add_argument("--ImvFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/AfterFiltering/mv.png"),
	help = "Location where the mean-variance image of the filtered dataset is stored")

parser.add_argument("--ImvNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Normalized/mv.png"),
	help = "Location where the mean-variance image of the normalized dataset is stored")

parser.add_argument("--ImvNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Normal/mv.png"),
	help = "Location where the mean-variance image of the dataset without transformed cells is stored")

parser.add_argument("--ImvWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/WithoutTissues/mv.png"),
	help = "Location where mean-variance images for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--ImvTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Top1000/mv.png"),
	help = "Location where the mean-variance image of the top 1000 dataset is stored")

parser.add_argument("--ImvTop88",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Top88/mv.png"),
	help = "Location where the mean-variance image of the top 88 dataset is stored")

parser.add_argument("--ImvCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/variants_chaperones/mv.png"),
	help = "Location where the mean-variance image of the variants and chaperones is stored")

parser.add_argument("--ImvTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/ByTissue/"),
	help = "Location where the mean-variance image of gene counts by tissue are stored")

parser.add_argument("--ImvRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Random/"),
	help = "Location where the mean-variance image of the random datasets is stored")

## PCA
parser.add_argument("--IpcaRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/BeforeFiltering/pca.png"),
	help = "Location where the pca image of the rawcounts is stored")

parser.add_argument("--IpcaFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/AfterFiltering/pca.png"),
	help = "Location where the pca image of the filtered non normalized counts is stored")

parser.add_argument("--IpcaNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/Normalized/pca.png"),
	help = "Location where the pca image of the filtered normalized counts is stored")

parser.add_argument("--IpcaNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/Normal/pca.png"),
	help = "Location where the pca image of the dataset without transformed cells is stored")

parser.add_argument("--IpcaWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/WithoutTissues/pca.png"),
	help = "Location where the pca image for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--IpcaTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/Top1000/pca.png"),
	help = "Location where the pca image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IpcaTop88",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/Top88/pca.png"),
	help = "Location where the pca image of the top 88 expressed gene counts is stored")

parser.add_argument("--IpcaCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/variants_chaperones/pca.png"),
	help = "Location where the pca image of the variant and chaperone gene counts is stored")

parser.add_argument("--IpcaTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/PCAbyTissue/"),
	help = "Location where the pca image of the gene counts by tissue are stored")

parser.add_argument("--IpcaRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/PCA/Random/pca.png"),
	help = "Location where the pca image of the randomly selected gene counts is stored")

## T-SNE
parser.add_argument("--ItsneRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/BeforeFiltering/T-Sne.png"),
	help = "Location where the T-Sne image of the rawcounts is stored")

parser.add_argument("--ItsneFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/AfterFiltering/T-Sne.png"),
	help = "Location where the T-Sne image of the filtered non normalized counts is stored")

parser.add_argument("--ItsneNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/Normalized/T-Sne.png"),
	help = "Location where the t-sne image of the filtered normalized counts is stored")

parser.add_argument("--ItsneNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/Normal/tsne.png"),
	help = "Location where the t-sne image of the dataset without transformed cells is stored")

parser.add_argument("--ItsneWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/WithoutTissues/tsne.png"),
	help = "Location where the t-sne image for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--ItsneTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/Top1000/T-Sne.png"),
	help = "Location where the t-sne image of the top 1000 expressed gene counts is stored")

parser.add_argument("--ItsneTop88",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/Top88/T-Sne.png"),
	help = "Location where the t-sne image of the top 88 expressed gene counts is stored")

parser.add_argument("--ItsneCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/variants_chaperones/T-Sne.png"),
	help = "Location where the t-sne image of the variant and chaperone gene counts is stored")

parser.add_argument("--ItsneTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/TSNEbyTissue/"),
	help = "Location where the t-sne image of the gene counts by tissue are stored")

parser.add_argument("--ItsneRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/T-Sne/Random/"),
	help = "Location where the T-Sne image of the randomly selected gene counts is stored")

## Metric
parser.add_argument("--distance",
	type = str,
	default = "euclidean",
	help = "Distance metric")

distance_metric = parser.parse_args().distance

## Clustermap
parser.add_argument("--IclstrNormS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples/" + distance_metric + "/Normalized/clustermap.png"),
	help = "Location where the clustermap of the normalized gene counts is stored")

parser.add_argument("--IclusterTopS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples/" + distance_metric + "/Top1000/clustermap.png"),
	help = "Location where the clustermap image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IclusterTop88S",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples/" + distance_metric + "/Top88/clustermap.png"),
	help = "Location where the clustermap image of the top 88 expressed gene counts is stored")

parser.add_argument("--IclusterCVs",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples/" + distance_metric + "/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap image of the variant and chaperone gene counts is stored")

parser.add_argument("--IclstrNonRcvS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples/" + distance_metric + "/nonReplicative/clustermap.png"),
	help = "Location where the clustermap of the histones chaperone and non replicative histone variant genes are stored")

parser.add_argument("--IclstrTissueS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples/" + distance_metric + "/ByTissue/"),
	help = "Location where the clustermap of the gene counts by tissue are stored")

parser.add_argument("--IclusterRandS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/" + distance_metric + "/Samples/Random/"),
	help = "Location where the clustermap image of the randomly selected gene counts is stored")

parser.add_argument("--IclstrNormG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Genes/" + distance_metric + "/Normalized/clustermap.png"),
	help = "Location where the clustermap of the normalized gene counts is stored")

parser.add_argument("--IclusterTopG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Genes/" + distance_metric + "/Top1000/clustermap.png"),
	help = "Location where the clustermap image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IclusterTop88G",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Genes/" + distance_metric + "/Top88/clustermap.png"),
	help = "Location where the clustermap image of the top 88 expressed gene counts is stored")

parser.add_argument("--IclusterCVG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Genes/" + distance_metric + "/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap image of the variant and chaperone gene counts is stored")

parser.add_argument("--IclstrNonRcvG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Genes/" + distance_metric + "/nonReplicative/clustermap.png"),
	help = "Location where the clustermap of the histones chaperone and non replicative histone variant genes are stored")

parser.add_argument("--IclstrTissueG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Genes/" + distance_metric + "/ByTissue/"),
	help = "Location where the clustermap of the gene counts by tissue are stored")

parser.add_argument("--IclusterRandG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Genes/" + distance_metric + "/Random/"),
	help = "Location where the clustermap image of the randomly selected gene counts is stored")

parser.add_argument("--IclstrNormSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples_Genes/" + distance_metric + "/Normalized/clustermap.png"),
	help = "Location where the clustermap of the normalized gene counts is stored")

parser.add_argument("--IclusterTopSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples_Genes/" + distance_metric + "/Top1000/clustermap.png"),
	help = "Location where the clustermap image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IclusterTop88SG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples_Genes/" + distance_metric + "/Top88/clustermap.png"),
	help = "Location where the clustermap image of the top 88 expressed gene counts is stored")

parser.add_argument("--IclusterCVsg",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples_Genes/" + distance_metric + "/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap image of the variant and chaperone gene counts is stored")

parser.add_argument("--IclstrNonRcvSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples_Genes/" + distance_metric + "/nonReplicative/clustermap.png"),
	help = "Location where the clustermap of the histones chaperone and non replicative histone variant genes are stored")

parser.add_argument("--IclusterTissueSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples_Genes/" + distance_metric + "/ByTissue/"),
	help = "Location where the clustermap image of the gene counts by tissue are stored")

parser.add_argument("--IclusterRandSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + choice.batch + "/Clustermap/Samples_Genes/" + distance_metric + "/Random/"),
	help = "Location where the clustermap image of the randomly selected gene counts is stored")

# Plotly
## General
parser.add_argument("--PgeneralRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/BeforeFiltering/"),
	help = "Location where QC html files for the raw dataset are stored")

parser.add_argument("--PgeneralFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/AfterFiltering/"),
	help = "Location where QC html files for the filtered dataset are stored")

parser.add_argument("--PgeneralNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/Normalized/"),
	help = "Location where QC html files for the normalized dataset are stored")

parser.add_argument("--PgeneralNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/Normal/"),
	help = "Location where QC images for the dataset without transformed cells are stored")

parser.add_argument("--PgeneralWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/WithoutTissues/"),
	help = "Location where QC html files for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--PgeneralTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/Top1000/"),
	help = "Location where QC html files for the top1000 expressed genes dataset are stored")

parser.add_argument("--PgeneralTop88",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/Top88/"),
	help = "Location where QC html files for the top 88 expressed genes dataset are stored")

parser.add_argument("--PgeneralCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/variants_chaperones/"),
	help = "Location where QC html files for the chaperone and variants dataset are stored")

parser.add_argument("--PgeneralNrCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/NonReplicativeHistones/"),
	help = "Location where QC html files for the chaperone and non replicative variants dataset are stored")

parser.add_argument("--PgeneralRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/General/Random/"),
	help = "Location where QC html files for the Random dataset are stored")

## Mean Variance
parser.add_argument("--PmvRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/BeforeFiltering/mv.html"),
	help = "Location where the mean-variance html file of the raw dataset is stored")

parser.add_argument("--PmvFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/AfterFiltering/mv.html"),
	help = "Location where the mean-variance html file of the filtered dataset is stored")

parser.add_argument("--PmvNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Normalized/mv.html"),
	help = "Location where the mean-variance html file of the normalized dataset is stored")

parser.add_argument("--PmvNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Normal/mv.html"),
	help = "Location where the mean-variance html file of the dataset without transformed cells is stored")

parser.add_argument("--PmvWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/WithoutTissues/mv.html"),
	help = "Location where mean-variance html file for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--PmvTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Top1000/mv.html"),
	help = "Location where the mean-variance html file of the top1000 dataset is stored")

parser.add_argument("--PmvTop88",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Top88/mv.html"),
	help = "Location where the mean-variance html file of the top 88 dataset is stored")

parser.add_argument("--PmvCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/variants_chaperones/mv.html"),
	help = "Location where the mean-variance html file of the variants and chaperones is stored")

parser.add_argument("--PmvTissue",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/ByTissue/"),
	help = "Location where the mean-variance html file of the gene counts by tissue are stored")

parser.add_argument("--PmvRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/MV_Plots/Random/"),
	help = "Location where the mean-variance html files of the random datasets is stored")

## PCA
parser.add_argument("--PpcaRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/BeforeFiltering/pca.html"),
	help = "Location where the pca html file of the raw counts is stored")

parser.add_argument("--PpcaFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/AfterFiltering/pca.html"),
	help = "Location where the pca html file of the filtered non normalized counts is stored")

parser.add_argument("--PpcaNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/Normalized/pca.html"),
	help = "Location where the pca html file of the filtered normalized counts is stored")

parser.add_argument("--PpcaNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/Normal/pca.html"),
	help = "Location where the pca html file of the dataset without transformed cells is stored")

parser.add_argument("--PpcaWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/WithoutTissues/pca.html"),
	help = "Location where the pca html file for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--PpcaTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/Top1000/pca.html"),
	help = "Location where the pca html file of the top 1000 expressed gene counts is stored")

parser.add_argument("--PpcaTop88",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/Top88/pca.html"),
	help = "Location where the pca html file of the top 88 expressed gene counts is stored")

parser.add_argument("--PpcaCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/variants_chaperones/pca.html"),
	help = "Location where the pca html file of the variant and chaperone gene counts is stored")

parser.add_argument("--PpcaTissue",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/PCAbyTissue/"),
	help = "Location where the pca html files of the gene counts by tissue are stored")

parser.add_argument("--PpcaRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/PCA/Random/"),
	help = "Location where the pca html file of the random generated files is stored")

## T-SNE
parser.add_argument("--PtsneRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/BeforeFiltering/T-Sne.html"),
	help = "Location where the t-sne html file of the raw counts is stored")

parser.add_argument("--PtsneFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/AfterFiltering/T-Sne.html"),
	help = "Location where the t-sne html file of the filtered non normalized counts is stored")

parser.add_argument("--PtsneNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/Normalized/T-Sne.html"),
	help = "Location where the t-sne html file of the filtered normalized counts is stored")

parser.add_argument("--PtsneNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/Normal/tsne.html"),
	help = "Location where the t-sne html file of the dataset without transformed cells is stored")

parser.add_argument("--PtsneWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/WithoutTissues/tsne.html"),
	help = "Location where the t-sne html file for the dataset without brain, blood, bone, pituitary, spleen and testis is stored")

parser.add_argument("--PtsneTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/Top1000/T-Sne.html"),
	help = "Location where the t-sne html fileof the top 1000 expressed gene counts is stored")

parser.add_argument("--PtsneTop88",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/Top88/T-Sne.html"),
	help = "Location where the t-sne html fileof the top 88 expressed gene counts is stored")

parser.add_argument("--PtsneCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/variants_chaperones/T-Sne.html"),
	help = "Location where the t-sne html file of the variant and chaperone gene counts is stored")

parser.add_argument("--PtsneTissue",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/TSNEbyTissue"),
	help = "Location where the t-sne html file of the gene counts by tissue are stored")

parser.add_argument("--PtsneRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + choice.batch + "/T-Sne/Random/"),
	help = "Location where the t-sne html file of the randomly selected gene counts is stored")

args = parser.parse_args()