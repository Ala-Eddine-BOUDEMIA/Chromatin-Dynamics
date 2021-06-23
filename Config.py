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

parser.add_argument("--bf",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/BeforeFiltering/PairedEndRounded.tsv"),
	help = "Location where the raw counts are stored")

parser.add_argument("--af",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Afteriltering/FilteredCPM10S36.tsv"),
	help = "Location where the filtered counts are stored")

parser.add_argument("--norm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Normalized/Normalized.tsv"),
	help = "Location where the normalized filtered counts are stored")

parser.add_argument("--top1000",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Top1000/Top1000Genes.tsv"),
	help = "Location where the counts of the top 1000 expressed genes are stored")

parser.add_argument("--rand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/Random/"),
	help = "Location where the random selected gene counts are stored")

parser.add_argument("--cv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Counts/variants_chaperones/variants_chaperones_counts.tsv"),
	help = "Location where the counts of the histones chaperone and histone variant genes are stored")

parser.add_argument("--list",
	type = Path,
	default = Path("Data/variants_chaperones/complete_list.csv"),
	help = "Location where the metadata file is stored")

parser.add_argument("--meta",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Metadata/" + choice.dataset + ".tsv"),
	help = "Location where the list of histone and chaperone genes is stored")

parser.add_argument("--corrNormG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Genes/Normalized/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the filtered normalized gene counts is stored")

parser.add_argument("--corrTopG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Genes/Top1000/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--corrCVg",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Genes/variants_chaperones/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the variant and chaperone gene counts is stored")

parser.add_argument("--corrRandG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Genes/Random/"),
	help = "Location where the correlation matrix of the random selected gene counts is stored")

parser.add_argument("--corrNormS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Sampels/Normalized/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the filtered normalized gene counts is stored")

parser.add_argument("--corrTopS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Sampels/Top1000/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the top 1000 expressed gene counts is stored")

parser.add_argument("--corrCVs",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Sampels/variants_chaperones/corr_matrix.tsv"),
	help = "Location where the correlation matrix of the variant and chaperone gene counts is stored")

parser.add_argument("--corrRandS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/CorrelationMatrix/Samples/Random/"),
	help = "Location where the correlation matrix of the random selected gene counts is stored")

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

parser.add_argument("--ImvNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Normalized/mv.png"),
	help = "Location where the mean-variance image of the normalized dataset is stored")

parser.add_argument("--ImvRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Random/mv.png"),
	help = "Location where the mean-variance image of the random datasets is stored")

parser.add_argument("--ImvTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Top1000/mv.png"),
	help = "Location where the mean-variance image of the top1000 dataset is stored")

parser.add_argument("--ImvCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/variants_chaperones/mv.png"),
	help = "Location where the mean-variance image of the variants and chaperones is stored")

parser.add_argument("--IpcaNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Normalized/pca.png"),
	help = "Location where the pca image of the filtered normalized counts is stored")

parser.add_argument("--IpcaTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/Top1000/pca.png"),
	help = "Location where the pca image of the top 1000 expressed gene counts is stored")

parser.add_argument("--IpcaCV",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/PCA/variants_chaperones/pca.png"),
	help = "Location where the pca image of the variant and chaperone gene counts is stored")

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

# Plotly
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

parser.add_argument("--PmvNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/Normalized/mv.png"),
	help = "Location where the mean-variance html file of the normalized dataset is stored")

parser.add_argument("--PmvRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/Random/mv.png"),
	help = "Location where the mean-variance html files of the random datasets is stored")

parser.add_argument("--PmvTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/MV_Plots/Top1000/mv.png"),
	help = "Location where the mean-variance html file of the top1000 dataset is stored")

parser.add_argument("--PmvCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/variants_chaperones/mv.png"),
	help = "Location where the mean-variance html file of the variants and chaperones is stored")

parser.add_argument("--PpcaNorm",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/MV_Plots/Normalized/pca.png"),
	help = "Location where the pca html file of the filtered normalized counts is stored")

parser.add_argument("--PpcaTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/Top1000/pca.png"),
	help = "Location where the pca html file of the top 1000 expressed gene counts is stored")

parser.add_argument("--PpcaCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/PCA/variants_chaperones/pca.png"),
	help = "Location where the pca html file of the variant and chaperone gene counts is stored")

parser.add_argument("--PtsneNorm",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/T-Sne/Normalized/T-Sne.png"),
	help = "Location where the t-sne html file of the filtered normalized counts is stored")

parser.add_argument("--PtsneTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/Top1000/T-Sne.png"),
	help = "Location where the t-sne html fileof the top 1000 expressed gene counts is stored")

parser.add_argument("--PtsneCV",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/T-Sne/variants_chaperones/T-Sne.png"),
	help = "Location where the t-sne html file of the variant and chaperone gene counts is stored")

args = parser.parse_args()