import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

# GTEx or TCGA
parser.add_argument("--dataset",
	type = str,
	default = "GTEx",
	help = "Dataset to use GTEx or TCGA")

# Data
choice = parser.parse_args()

## Counts
parser.add_argument("--bf",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/BeforeFiltering/PairedEndRounded.tsv"),
	help = "Location where the raw counts are stored")

parser.add_argument("--af",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/AfterFiltering/FilteredCPM10S36.tsv"),
	help = "Location where the filtered counts are stored")

parser.add_argument("--norm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Normalized/Normalized.tsv"),
	help = "Location where the normalized filtered counts are stored")

parser.add_argument("--top1000",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Top1000/Top1000Genes.tsv"),
	help = "Location where the counts of the top 1000 expressed genes are stored")

parser.add_argument("--top76",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Top76/Top76Genes.tsv"),
	help = "Location where the counts of the top 76 expressed genes are stored")

parser.add_argument("--rand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Random/"),
	help = "Location where the random selected gene counts are stored")

parser.add_argument("--cv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/variants_chaperones/variants_chaperones_counts.tsv"),
	help = "Location where the counts of the histones chaperone and histone variant genes are stored")

parser.add_argument("--nonRcv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/variants_chaperones/nonRvariants_chaperones_counts.tsv"),
	help = "Location where the counts of the histones chaperone and non replicative histone variant genes are stored")

## Metadata
parser.add_argument("--list",
	type = Path,
	default = Path("Data/variants_chaperones/complete_list.csv"),
	help = "Location where the list of histone and chaperone genes is stored")

parser.add_argument("--nonReplicative",
	type = Path,
	default = Path("Data/variants_chaperones/withoutR.txt"),
	help = "Location where the list of non replicative histone genes and chaperone genes is stored")

parser.add_argument("--meta",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Metadata/" + choice.dataset + ".tsv"),
	help = "Location where the metadata file is stored")

## Spearman or Pearson
parser.add_argument("--corrMethod",
	type = str,
	default = "spearman",
	help = "The method by which the correlation matrix will be computed")

corr_method = parser.parse_args()
corr_method = corr_method.corrMethod

## Correlation
parser.add_argument("--corrNormG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Genes/Normalized/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the filtered normalized gene counts is stored")

parser.add_argument("--corrTopG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Genes/Top1000/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--corr76G",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Genes/Top76/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the top 76 expressed genes are stored")

parser.add_argument("--corrCVg",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Genes/variants_chaperones/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the variant and chaperone gene counts is stored")

parser.add_argument("--corrNonRcvG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Genes/nonReplicative/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the histones chaperone and non replicative histone variant genes are stored")

parser.add_argument("--corrRandG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Genes/Random/"),
	help = "Location where the correlation matrix of the random selected gene counts is stored")

parser.add_argument("--corrNormS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Samples/Normalized/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the filtered normalized gene counts is stored")

parser.add_argument("--corrTopS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Samples/Top1000/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--corr76S",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Samples/Top76/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the top 76 expressed genes are stored")

parser.add_argument("--corrCVs",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Samples/variants_chaperones/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the variant and chaperone gene counts is stored")

parser.add_argument("--corrNonRcvS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Samples/nonReplicative/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the histones chaperone and non replicative histone variant genes are stored")

parser.add_argument("--corrRandS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/" + corr_method + "/Samples/Random/"),
	help = "Location where the correlation matrix of the random selected gene counts is stored")

## PCA
parser.add_argument("--pcaRaw",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/PCA/BeforeFiltering/pca.tsv"),
	help = "Location where the pca matrix of the raw counts is stored")

parser.add_argument("--pcaFiltered",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/PCA/AfterFiltering/pca.tsv"),
	help = "Location where the pca matrix of the filtered non normalized counts is stored")

parser.add_argument("--pcaRand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/PCA/Random/"),
	help = "Location where the pca matrix of the randomly selected geen counts is stored")

parser.add_argument("--pcaNorm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/PCA/Normalized/pca.tsv"),
	help = "Location where the pca matrix of the filtered normalized counts is stored")

parser.add_argument("--pcaTop",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/PCA/Top1000/pca.tsv"),
	help = "Location where the pca matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--pcaCV",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/PCA/variants_chaperones/pca.tsv"),
	help = "Location where the pca matrix of the variant and chaperone gene counts is stored")

## T-SNE
parser.add_argument("--tsneRaw",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/T-Sne/BeforeFiltering/T-Sne.tsv"),
	help = "Location where the T-Sne matrix of the raw counts is stored")

parser.add_argument("--tsneFiltered",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/T-Sne/AfterFiltering/T-Sne.tsv"),
	help = "Location where the T-Sne matrix of the filtered non normalized counts is stored")

parser.add_argument("--tsneRand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/T-Sne/Random/"),
	help = "Location where the T-Sne matrix of the randomly selected geen counts is stored")

parser.add_argument("--tsneNorm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/T-Sne/Normalized/T-Sne.tsv"),
	help = "Location where the t-sne matrix of the filtered normalized counts is stored")

parser.add_argument("--tsneTop",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/T-Sne/Top1000/T-Sne.tsv"),
	help = "Location where the t-sne matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--tsneCV",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/T-Sne/variants_chaperones/T-Sne.tsv"),
	help = "Location where the t-sne matrix of the variant and chaperone gene counts is stored")

# Images
## General
parser.add_argument("--IgeneralRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/General/BeforeFiltering/"),
	help = "Location where QC images for the raw dataset are stored")

parser.add_argument("--IgeneralFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/General/AfterFiltering/"),
	help = "Location where QC images for the filtered dataset are stored")

parser.add_argument("--IgeneralNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/General/Normalized/"),
	help = "Location where QC images for the normalized dataset are stored")

parser.add_argument("--IgeneralRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/General/Random/"),
	help = "Location where QC images for the Random dataset are stored")

parser.add_argument("--IgeneralTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/General/Top1000/"),
	help = "Location where QC images for the top1000 expressed genes dataset are stored")

parser.add_argument("--IgeneralCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/General/variants_chaperones/"),
	help = "Location where QC images for the chaperone and variants dataset are stored")

parser.add_argument("--IgeneralNrCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/General/NonReplicativeHistones/"),
	help = "Location where QC images for the chaperone and non replicative variants dataset are stored")

## Mean Variance
parser.add_argument("--ImvRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/BeforeFiltering/mv.png"),
	help = "Location where the mean-variance image of the raw dataset is stored")

parser.add_argument("--ImvFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/AfterFiltering/mv.png"),
	help = "Location where the mean-variance image of the filtered dataset is stored")

parser.add_argument("--ImvNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Normalized/mv.png"),
	help = "Location where the mean-variance image of the normalized dataset is stored")

parser.add_argument("--ImvRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Random/"),
	help = "Location where the mean-variance image of the random datasets is stored")

parser.add_argument("--ImvTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Top1000/mv.png"),
	help = "Location where the mean-variance image of the top1000 dataset is stored")

parser.add_argument("--ImvCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/variants_chaperones/mv.png"),
	help = "Location where the mean-variance image of the variants and chaperones is stored")

## PCA
parser.add_argument("--IpcaRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/BeforeFiltering/pca.png"),
	help = "Location where the pca image of the rawcounts is stored")

parser.add_argument("--IpcaFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/AfterFiltering/pca.png"),
	help = "Location where the pca image of the filtered non normalized counts is stored")

parser.add_argument("--IpcaRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/Random/"),
	help = "Location where the pca image of the randomly selected gene counts is stored")

parser.add_argument("--IpcaNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/Normalized/pca.png"),
	help = "Location where the pca image of the filtered normalized counts is stored")

parser.add_argument("--IpcaTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/Top1000/pca.png"),
	help = "Location where the pca image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IpcaCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/variants_chaperones/pca.png"),
	help = "Location where the pca image of the variant and chaperone gene counts is stored")

## T-SNE
parser.add_argument("--ItsneRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/T-Sne/BeforeFiltering/T-Sne.png"),
	help = "Location where the T-Sne image of the rawcounts is stored")

parser.add_argument("--ItsneFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/T-Sne/AfterFiltering/T-Sne.png"),
	help = "Location where the T-Sne image of the filtered non normalized counts is stored")

parser.add_argument("--ItsneRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/T-Sne/Random/"),
	help = "Location where the T-Sne image of the randomly selected gene counts is stored")

parser.add_argument("--ItsneNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/T-Sne/Normalized/T-Sne.png"),
	help = "Location where the t-sne image of the filtered normalized counts is stored")

parser.add_argument("--ItsneTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/T-Sne/Top1000/T-Sne.png"),
	help = "Location where the t-sne image of the top 1000 expressed gene counts is stored")

parser.add_argument("--ItsneCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/T-Sne/variants_chaperones/T-Sne.png"),
	help = "Location where the t-sne image of the variant and chaperone gene counts is stored")

## Clustermap
parser.add_argument("--IclusterRandS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Samples/Random/"),
	help = "Location where the clustermap image of the randomly selected gene counts is stored")

parser.add_argument("--IclusterNormS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "Samples/Normalized/clustermap.png"),
	help = "Location where the clustermap image of the filtered normalized counts is stored")

parser.add_argument("--IclusterTopS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Samples/Top1000/clustermap.png"),
	help = "Location where the clustermap image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IclusterTop76S",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Samples/Top76/clustermap.png"),
	help = "Location where the clustermap image of the top 76 expressed gene counts is stored")

parser.add_argument("--IclusterCVs",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Samples/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap image of the variant and chaperone gene counts is stored")

parser.add_argument("--IclstrNonRcvS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Samples/nonReplicative/clustermap.png"),
	help = "Location where the clustermap of the histones chaperone and non replicative histone variant genes are stored")

parser.add_argument("--IclusterRandG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Genes/Random/"),
	help = "Location where the clustermap image of the randomly selected gene counts is stored")

parser.add_argument("--IclusterNormG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Genes/Normalized/clustermap.png"),
	help = "Location where the clustermap image of the filtered normalized counts is stored")

parser.add_argument("--IclusterTopG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Genes/Top1000/clustermap.png"),
	help = "Location where the clustermap image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IclusterTop76G",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Genes/Top76/clustermap.png"),
	help = "Location where the clustermap image of the top 76 expressed gene counts is stored")

parser.add_argument("--IclstrNonRcvG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Genes/nonReplicative/clustermap.png"),
	help = "Location where the clustermap of the histones chaperone and non replicative histone variant genes are stored")

parser.add_argument("--IclusterCVG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/" + corr_method + "/Genes/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap image of the variant and chaperone gene counts is stored")

parser.add_argument("--IclusterRandSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/Random/"),
	help = "Location where the clustermap image of the randomly selected gene counts is stored")

parser.add_argument("--IclusterNormSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/Normalized/clustermap.png"),
	help = "Location where the clustermap image of the filtered normalized counts is stored")

parser.add_argument("--IclusterTopSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/Top1000/clustermap.png"),
	help = "Location where the clustermap image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IclusterTop76SG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/Top76/clustermap.png"),
	help = "Location where the clustermap image of the top 76 expressed gene counts is stored")

parser.add_argument("--IclusterCVsg",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap image of the variant and chaperone gene counts is stored")

parser.add_argument("--IclstrNonRcvSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/nonReplicative/clustermap.png"),
	help = "Location where the clustermap of the histones chaperone and non replicative histone variant genes are stored")

# Plotly
## General
parser.add_argument("--PgeneralRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/General/BeforeFiltering/"),
	help = "Location where QC html files for the raw dataset are stored")

parser.add_argument("--PgeneralFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/General/AfterFiltering/"),
	help = "Location where QC html files for the filtered dataset are stored")

parser.add_argument("--PgeneralNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/General/Normalized/"),
	help = "Location where QC html files for the normalized dataset are stored")

parser.add_argument("--PgeneralRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/General/Random/"),
	help = "Location where QC html files for the Random dataset are stored")

parser.add_argument("--PgeneralTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/General/Top1000/"),
	help = "Location where QC html files for the top1000 expressed genes dataset are stored")

parser.add_argument("--PgeneralCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/General/variants_chaperones/"),
	help = "Location where QC html files for the chaperone and variants dataset are stored")

parser.add_argument("--PgeneralNrCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/General/NonReplicativeHistones/"),
	help = "Location where QC html files for the chaperone and non replicative variants dataset are stored")

## Mean Variance
parser.add_argument("--PmvRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/BeforeFiltering/mv.html"),
	help = "Location where the mean-variance html file of the raw dataset is stored")

parser.add_argument("--PmvFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/AfterFiltering/mv.html"),
	help = "Location where the mean-variance html file of the filtered dataset is stored")

parser.add_argument("--PmvNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/Normalized/mv.html"),
	help = "Location where the mean-variance html file of the normalized dataset is stored")

parser.add_argument("--PmvRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/Random/"),
	help = "Location where the mean-variance html files of the random datasets is stored")

parser.add_argument("--PmvTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/Top1000/mv.html"),
	help = "Location where the mean-variance html file of the top1000 dataset is stored")

parser.add_argument("--PmvCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/variants_chaperones/mv.html"),
	help = "Location where the mean-variance html file of the variants and chaperones is stored")

## PCA
parser.add_argument("--PpcaRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/BeforeFiltering/pca.html"),
	help = "Location where the pca html file of the raw counts is stored")

parser.add_argument("--PpcaFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/AfterFiltering/pca.html"),
	help = "Location where the pca html file of the filtered non normalized counts is stored")

parser.add_argument("--PpcaRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/Random/"),
	help = "Location where the pca html file of the random generated files is stored")

parser.add_argument("--PpcaNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/Normalized/pca.html"),
	help = "Location where the pca html file of the filtered normalized counts is stored")

parser.add_argument("--PpcaTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/Top1000/pca.html"),
	help = "Location where the pca html file of the top 1000 expressed gene counts is stored")

parser.add_argument("--PpcaCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/variants_chaperones/pca.html"),
	help = "Location where the pca html file of the variant and chaperone gene counts is stored")

## T-SNE
parser.add_argument("--PtsneRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/BeforeFiltering/T-Sne.html"),
	help = "Location where the t-sne html file of the raw counts is stored")

parser.add_argument("--PtsneFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/AfterFiltering/T-Sne.html"),
	help = "Location where the t-sne html file of the filtered non normalized counts is stored")

parser.add_argument("--PtsneRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/Random/"),
	help = "Location where the t-sne html file of the randomly selected gene counts is stored")

parser.add_argument("--PtsneNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/Normalized/T-Sne.html"),
	help = "Location where the t-sne html file of the filtered normalized counts is stored")

parser.add_argument("--PtsneTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/Top1000/T-Sne.html"),
	help = "Location where the t-sne html fileof the top 1000 expressed gene counts is stored")

parser.add_argument("--PtsneCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/variants_chaperones/T-Sne.html"),
	help = "Location where the t-sne html file of the variant and chaperone gene counts is stored")

## Clustermap
parser.add_argument("--PclusterRandS",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples/Random/"),
	help = "Location where the clustermap html file of the randomly selected gene counts is stored")

parser.add_argument("--PclusterNormS",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples/Normalized/clustermap.png"),
	help = "Location where the clustermap html file of the filtered normalized counts is stored")

parser.add_argument("--PclusterTopS",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples/Top1000/clustermap.png"),
	help = "Location where the clustermap html file of the top 1000 expressed gene counts is stored")

parser.add_argument("--PclusterCVs",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap html file of the variant and chaperone gene counts is stored")

parser.add_argument("--PclusterRandG",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Genes/Random/"),
	help = "Location where the clustermap html file of the randomly selected gene counts is stored")

parser.add_argument("--PclusterNormG",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Genes/Normalized/clustermap.png"),
	help = "Location where the clustermap html file of the filtered normalized counts is stored")

parser.add_argument("--PclusterTopG",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Genes/Top1000/clustermap.png"),
	help = "Location where the clustermap html file of the top 1000 expressed gene counts is stored")

parser.add_argument("--PclusterCVG",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Genes/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap html file of the variant and chaperone gene counts is stored")

parser.add_argument("--PclusterRandSG",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/Random/"),
	help = "Location where the clustermap html file of the randomly selected gene counts is stored")

parser.add_argument("--PclusterNormSG",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/Normalized/clustermap.png"),
	help = "Location where the clustermap html file of the filtered normalized counts is stored")

parser.add_argument("--PclusterTopSG",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/Top1000/clustermap.png"),
	help = "Location where the clustermap html file of the top 1000 expressed gene counts is stored")

parser.add_argument("--PclusterCVsg",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/Clustermap/Samples_Genes/variants_chaperones/clustermap.png"),
	help = "Location where the clustermap html file of the variant and chaperone gene counts is stored")

args = parser.parse_args()