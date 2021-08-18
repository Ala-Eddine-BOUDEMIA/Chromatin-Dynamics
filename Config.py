import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

# GTEx or TCGA or GTEx_TCGA
parser.add_argument("--dataset",
	type = str,
	default = "TCGA",
	help = "Dataset to use: GTEx or TCGA or GTEx_TCGA")

# Method used to normalize (CPM or TMM)
parser.add_argument("--normMethod",
	type = str,
	default = "CPM",
	help = "Method to normalize")

# Data to be used for the analysis
# Either the full normalized dataset
# Or the histone chaperones and histone variants dataset
# Possible arguments: Normalized - variants_chaperones
parser.add_argument("--which",
	type = str,
	default = "variants_chaperones",
	help = "Dataset to use: normalized or \
			histone chaperones and histone variants")

# Data
choice = parser.parse_args()

# tissues or cancer
if choice.dataset == "GTEx" or choice.dataset == "GTEx_TCGA":
	var1 = "smtsd"
	var2 = "smts"
else:
	var1 = "gdc_cases.tissue_source_site.project"
	var2 = "gdc_cases.tissue_source_site.project"

parser.add_argument("--smtsd",
	type = str,
	default = var1,
	help = "-")

parser.add_argument("--smts",
	type = str,
	default = var2,
	help = "-")

# bigwig or run
if choice.dataset == "GTEx" or choice.dataset == "GTEx_TCGA":
	var = "run"
else:
	var = "bigwig_file"

parser.add_argument("--id",
	type = str,
	default = var,
	help = "-")

## Counts
parser.add_argument("--bf",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset \
			+ "/BeforeFiltering/PairedEndRounded.tsv"),
	help = "Raw counts")

parser.add_argument("--af",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/AfterFiltering/FilteredCPM5S18.tsv"),
	help = "Filtered counts, you can change the file name to one of these: \
			Filtered + CPM5S18 or CPM10S18 or CPM10S36 + .tsv")

parser.add_argument("--norm",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/Normalized/Full/counts.tsv"),
	help = "Normalized counts")

parser.add_argument("--full",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/" + choice.which + "/Full/counts.tsv"),
	help = choice.which + " complete counts")

parser.add_argument("--onlyNormal",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/" + choice.which + "/Normal/counts.tsv"),
	help = "Counts missing the samples coming from transformed cells")

parser.add_argument("--WoTissues",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/" + choice.which + "/WithoutTissues/counts.tsv"),
	help = "Counts missing samples coming from:\
			Blood, Brain, Bone, Pituitary and Spleen tissues")

parser.add_argument("--top1000",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/Top1000/Top1000Genes.tsv"),
	help = "Counts limited to the Top 1000 expressed genes")

parser.add_argument("--top100",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/Top100/Top100Genes.tsv"),
	help = "Counts limited to the Top 100 (88-127) expressed genes")

parser.add_argument("--nonRcv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/variants_chaperones/NonReplicative/counts.tsv"),
	help = "Counts limited to the histone chaperones \
			and histone non-replicative variants genes")

parser.add_argument("--tissue",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/" + choice.which + "/CountsByTissue/"),
	help = "Counts by tissue type")

parser.add_argument("--rand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Counts/Random/"),
	help = "Randomly selected counts")

## Metadata
parser.add_argument("--meta",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/Metadata/" \
			+ choice.dataset + ".tsv"),
	help = "Metadata file")

parser.add_argument("--list",
	type = Path,
	default = Path("Data/variants_chaperones/complete_list.csv"),
	help = "List of the names and IDs of histone chaperones \
			and histone variants genes")

parser.add_argument("--nonReplicative",
	type = Path,
	default = Path("Data/variants_chaperones/withoutR.txt"),
	help = "List of the names and IDs of histone chaperones \
			and histone non-replicative variants genes")

## Correlation
parser.add_argument("--corrFullG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/" + choice.which \
			+ "/Full/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrNormalG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/" + choice.which \
			+ "/Normal/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrWoTissuesG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/" + choice.which \
			+ "/WithoutTissues/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrTopG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/Top1000/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrTop100G",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/Top100/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrNRcvG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/variants_chaperones/NonReplicative/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrTissueG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/" + choice.which \
			+ "/CorrByTissue/"),
	help = "Correlation matrices")

parser.add_argument("--corrRandG",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Genes/Random/"),
	help = "Correlation matrices")

parser.add_argument("--corrFullS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Samples/" + choice.which \
			+ "/Full/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrNormalS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Samples/" + choice.which \
			+ "/Normal/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrWoTissuesS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Samples/" + choice.which \
			+ "/WithoutTissues/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrTopS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
		+ "/CorrelationMatrix/Samples/Top1000/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrTop100S",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Samples/Top100/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrNRcvS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Samples/variants_chaperones/NonReplicative/corr_matrix.tsv"),
	help = "Correlation matrix")

parser.add_argument("--corrTissueS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Samples/" + choice.which + "/CorrByTissue/"),
	help = "Correlation matrices")

parser.add_argument("--corrRandS",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/CorrelationMatrix/Samples/Random/"),
	help = "Correlation matrices")

## PCA
parser.add_argument("--pcaRaw",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/BeforeFiltering/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaFiltered",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/AfterFiltering/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaFull",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/Full/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaNormal",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/Normal/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaWoTissues",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/WithoutTissues/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaTop",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/Top1000/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaTop100",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/Top100/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaNRcv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/variants_chaperones/NonReplicative/pca.tsv"),
	help = "PCA file")

parser.add_argument("--pcaTissue",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/PCAbyTissue/"),
	help = "PCA files")

parser.add_argument("--pcaRand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/Random/"),
	help = "PCA file")

## T-SNE
parser.add_argument("--tsneRaw",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/BeforeFiltering/T-Sne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneFiltered",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/AfterFiltering/T-Sne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneFull",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/Full/tsne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneNormal",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/Normal/tsne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneWoTissues",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/WithoutTissues/tsne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneTop",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/Top1000/T-Sne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneTop100",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/Top100/T-Sne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneNRcv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/variants_chaperones/NonReplicative/T-Sne.tsv"),
	help = "T-SNE file")

parser.add_argument("--tsneTissue",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/TSNEbyTissue/"),
	help = "T-SNE file")

parser.add_argument("--tsneRand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/Random/"),
	help = "T-SNE file")

## UMAP
parser.add_argument("--umapRaw",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/BeforeFiltering/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapFiltered",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/AfterFiltering/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapFull",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/Full/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapNormal",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/Normal/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapWoTissues",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/WithoutTissues/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapTop",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/Top1000/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapTop100",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/Top100/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapNRcv",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/variants_chaperones/NonReplicative/umap.tsv"),
	help = "UMAP file")

parser.add_argument("--umapTissue",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/umapbyTissue/"),
	help = "UMAP file")

parser.add_argument("--umapRand",
	type = Path,
	default = Path("Data/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/Random/"),
	help = "UMAP file")

# Images
## General
parser.add_argument("--IgeneralRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			 + "/General/BeforeFiltering/"),
	help = "QC images")

parser.add_argument("--IgeneralFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/AfterFiltering/"),
	help = "QC images")

parser.add_argument("--IgeneralFull",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/" + choice.which + "/Full/"),
	help = "QC images")

parser.add_argument("--IgeneralNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/" + choice.which + "/Normal/"),
	help = "QC images")

parser.add_argument("--IgeneralWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/" + choice.which + "/WithoutTissues/"),
	help = "QC images")

parser.add_argument("--IgeneralTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/Top1000/"),
	help = "QC images")

parser.add_argument("--IgeneralTop100",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/Top100/"),
	help = "QC images")

parser.add_argument("--IgeneralNRcv",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/variants_chaperones/NonReplicative/"),
	help = "QC images")

parser.add_argument("--IgeneralRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/General/Random/"),
	help = "QC images")

## Mean Variance
parser.add_argument("--ImvRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/BeforeFiltering/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/AfterFiltering/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvFull",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/" + choice.which + "/Full/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/" + choice.which + "/Normal/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/" + choice.which + "/WithoutTissues/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/Top1000/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvTop100",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/Top100/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvNRcv",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/variants_chaperones/NonReplicative/mv.png"),
	help = "Mean-variance image")

parser.add_argument("--ImvTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/" + choice.which + "/ByTissue/"),
	help = "Mean-variance images")

parser.add_argument("--ImvRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/MV_Plots/Random/"),
	help = "Mean-variance images")

## Z scores
parser.add_argument("--IzscoreRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/BeforeFiltering/Z_scores.png"),
	help = "Z_scores image")

parser.add_argument("--IzscoreFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/AfterFiltering/Z_scores.png"),
	help = "Z_scores image")

parser.add_argument("--IzscoreFull",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/" + choice.which + "/Full/Z_scores.png"),
	help = "Z_scores image")

parser.add_argument("--IzscoreNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/" + choice.which + "/Normal/Z_scores.png"),
	help = "Z_scores image")

parser.add_argument("--IzscoreWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/" + choice.which + "/WithoutTissues/Z_scores.png"),
	help = "Z_scores images")

parser.add_argument("--IzscoreTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/Top1000/Z_scores.png"),
	help = "Z_scores image")

parser.add_argument("--IzscoreTop100",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/Top100/Z_scores.png"),
	help = "Z_scores image")

parser.add_argument("--IzscoreNRcv",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/variants_chaperones/NonReplicative/Z_scores.png"),
	help = "Z_scores image")

parser.add_argument("--IzscoreTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/" + choice.which + "/ByTissue/"),
	help = "Z_scores images")

parser.add_argument("--IzscoreRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Z_scores/Random/"),
	help = "Z_scores images")

## PCA
parser.add_argument("--IpcaRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/BeforeFiltering/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/AfterFiltering/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaFull",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/Full/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/Normal/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/WithoutTissues/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/Top1000/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaTop100",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/Top100/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaNRcv",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/variants_chaperones/NonReplicative/pca.png"),
	help = "PCA image")

parser.add_argument("--IpcaTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/" + choice.which + "/PCAbyTissue/"),
	help = "PCA images")

parser.add_argument("--IpcaRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/PCA/Random/pca.png"),
	help = "PCA images")

## T-SNE
parser.add_argument("--ItsneRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/BeforeFiltering/T-Sne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/AfterFiltering/T-Sne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneFull",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/Full/T-Sne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/Normal/tsne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/WithoutTissues/tsne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/Top1000/T-Sne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneTop100",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/Top100/T-Sne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneNRcv",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/variants_chaperones/NonReplicative/T-Sne.png"),
	help = "T-SNE image")

parser.add_argument("--ItsneTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/" + choice.which + "/TSNEbyTissue/"),
	help = "T-SNE images")

parser.add_argument("--ItsneRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/T-Sne/Random/"),
	help = "T-SNE images")

## UMAP
parser.add_argument("--IumapRaw",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/BeforeFiltering/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapFiltered",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/AfterFiltering/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapFull",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/Full/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapNormal",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/Normal/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapWoTissues",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/WithoutTissues/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapTop",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/Top1000/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapTop100",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/Top100/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapNRcv",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/variants_chaperones/NonReplicative/umap.png"),
	help = "UMAP image")

parser.add_argument("--IumapTissue",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/" + choice.which + "/umapbyTissue/"),
	help = "UMAP images")

parser.add_argument("--IumapRand",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/UMAP/Random/"),
	help = "UMAP images")

## Metric
parser.add_argument("--distance",
	type = str,
	default = "euclidean",
	help = "Distance metric")

distance_metric = parser.parse_args().distance

## Clustermap
parser.add_argument("--IclstrFullS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric + "/" + choice.which \
			+ "/Full/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrNormalS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric + "/" +choice.which \
			+ "/Normal/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrWoTissueS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric + "/" + choice.which \
			+ "/WithoutTissues/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTopS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric \
			+ "/Top1000/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTop100S",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric \
			+ "/Top100/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrNRcvS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric \
			+ "/variants_chaperones/NonReplicative/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTissueS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric + "/" + choice.which \
			+ "/ByTissue/"),
	help = "Clustermaps")

parser.add_argument("--IclstrRandS",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples/" + distance_metric + "/Random/"),
	help = "Clustermaps")

parser.add_argument("--IclstrFullG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric + "/" + choice.which \
			+ "/Full/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrNormalG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric + "/" + choice.which \
			+ "/Normal/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrWoTissueG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric + "/" + choice.which \
			+ "/WithoutTissues/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTopG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric \
			+ "/Top1000/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTop100G",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric + "/Top100/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrNRcvG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric \
			+ "/variants_chaperones/NonReplicative/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTissueG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric + "/" + choice.which \
			+ "/ByTissue/"),
	help = "Clustermaps")

parser.add_argument("--IclstrRandG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Genes/" + distance_metric + "/Random/"),
	help = "Clustermaps")

parser.add_argument("--IclstrFullSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric + "/" + choice.which \
			+ "/Full/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrNormalSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric + "/" + choice.which \
			+ "/Normal/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrWoTissueSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric + "/" + choice.which \
			+ "/WithoutTissues/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTopSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric \
			+ "/Top1000/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTop100SG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric + \
			"/Top100/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrNRcvSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric \
			+ "/variants_chaperones/NonReplicative/clustermap.png"),
	help = "Clustermap")

parser.add_argument("--IclstrTissueSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric + "/" + choice.which \
			+ "/ByTissue/"),
	help = "Clustermaps")

parser.add_argument("--IclstrRandSG",
	type = Path,
	default = Path("Images/").joinpath(choice.dataset + "/" + choice.normMethod \
			+ "/Clustermap/Samples_Genes/" + distance_metric + "/Random/"),
	help = "Clustermaps")

# Plotly
## General
parser.add_argument("--PgeneralRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/BeforeFiltering/"),
	help = "QC html files")

parser.add_argument("--PgeneralFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/AfterFiltering/"),
	help = "QC html files")

parser.add_argument("--PgeneralFull",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/" + choice.which + "/Full/"),
	help = "QC html files")

parser.add_argument("--PgeneralNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/" + choice.which + "/Normal/"),
	help = "QC html files")

parser.add_argument("--PgeneralWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/" + choice.which \
			+ "/WithoutTissues/"),
	help = "QC html files")

parser.add_argument("--PgeneralTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/Top1000/"),
	help = "QC html files")

parser.add_argument("--PgeneralTop100",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/Top100/"),
	help = "QC html files")

parser.add_argument("--PgeneralNRcv",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/variants_chaperones/" \
			+ "NonReplicative/"),
	help = "QC html files")

parser.add_argument("--PgeneralRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/General/Random/"),
	help = "QC html files")

## Mean Variance
parser.add_argument("--PmvRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/BeforeFiltering/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/AfterFiltering/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvFull",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/" + choice.which \
			+ "/Full/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/" + choice.which \
			+ "/Normal/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/" + choice.which \
			+ "/WithoutTissues/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/Top1000/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvTop100",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/Top100/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvNRcv",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/variants_chaperones/" \
			+ "NonReplicative/mv.html"),
	help = "Mean-variance html file")

parser.add_argument("--PmvTissue",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/" + choice.which \
			+ "/ByTissue/"),
	help = "Mean-variance html files")

parser.add_argument("--PmvRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset \
			+ "/" + choice.normMethod + "/MV_Plots/Random/"),
	help = "Mean-variance html files")

# Z scores
parser.add_argument("--PzscoreRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/BeforeFiltering/Z_scores.html"),
	help = "Z_scores html file")

parser.add_argument("--PzscoreFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/AfterFiltering/Z_scores.html"),
	help = "Z_scores html file")

parser.add_argument("--PzscoreFull",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			 + choice.normMethod + "/Z_scores/" + choice.which \
			 + "/Full/Z_scores.html"),
	help = "Z_scores html file")

parser.add_argument("--PzscoreNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/" + choice.which \
			+ "/Normal/Z_scores.html"),
	help = "Z_scores html file")

parser.add_argument("--PzscoreWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/" + choice.which \
			+ "/WithoutTissues/Z_scores.html"),
	help = "Z_scores images")

parser.add_argument("--PzscoreTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/Top1000/Z_scores.html"),
	help = "Z_scores html file")

parser.add_argument("--PzscoreTop100",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/Top100/Z_scores.html"),
	help = "Z_scores html file")

parser.add_argument("--PzscoreNRcv",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/variants_chaperones/" \
			+ "NonReplicative/Z_scores.html"),
	help = "Z_scores image")

parser.add_argument("--PzscoreTissue",
	type = Path,
	default = Path("Plotly_HTML_Files").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/" + choice.which + "/ByTissue/"),
	help = "Z_scores html file")

parser.add_argument("--PzscoreRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/Z_scores/Random/"),
	help = "Z_scores html file")

## PCA
parser.add_argument("--PpcaRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/BeforeFiltering/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/AfterFiltering/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaFull",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/" + choice.which + "/Full/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			 + choice.normMethod + "/PCA/" + choice.which + "/Normal/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/" + choice.which + "/WithoutTissues/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/Top1000/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaTop100",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/Top100/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaNRcv",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/variants_chaperones/" \
			+ "NonReplicative/pca.html"),
	help = "PCA html file")

parser.add_argument("--PpcaTissue",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/" + choice.which + "/PCAbyTissue/"),
	help = "PCA html files")

parser.add_argument("--PpcaRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/PCA/Random/"),
	help = "PCA html files")

## T-SNE
parser.add_argument("--PtsneRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/BeforeFiltering/T-Sne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/AfterFiltering/T-Sne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneFull",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/" + choice.which + "/Full/T-Sne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/" + choice.which + "/Normal/tsne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/" + choice.which \
			+ "/WithoutTissues/tsne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/Top1000/T-Sne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneTop100",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/Top100/T-Sne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneNRcv",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/variants_chaperones/" \
			+ "NonReplicative/T-Sne.html"),
	help = "T-SNE html file")

parser.add_argument("--PtsneTissue",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/" + choice.which + "/TSNEbyTissue"),
	help = "T-SNE html file")

parser.add_argument("--PtsneRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/T-Sne/Random/"),
	help = "T-SNE html file")

## UMAP
parser.add_argument("--PumapRaw",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/BeforeFiltering/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapFiltered",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/AfterFiltering/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapFull",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/" + choice.which + "/Full/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapNormal",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/" + choice.which + "/Normal/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapWoTissues",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/" + choice.which \
			+ "/WithoutTissues/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapTop",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/Top1000/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapTop100",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/Top100/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapNRcv",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/variants_chaperones/" \
			+ "NonReplicative/umap.html"),
	help = "UMAP html file")

parser.add_argument("--PumapTissue",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/" + choice.which + "/umapbyTissue"),
	help = "UMAP html file")

parser.add_argument("--PumapRand",
	type = Path,
	default = Path("Plotly_HTML_Files/").joinpath(choice.dataset + "/" \
			+ choice.normMethod + "/UMAP/Random/"),
	help = "UMAP html file")

args = parser.parse_args()

### Get counts
counts = [args.full, args.nonRcv]
""" [args.full, args.onlyNormal, args.WoTissues, 
	args.top1000, args.top100, args.nonRcv]"""

### Get genes' correlations
g_corr = [args.corrFullG, args.corrNRcvG] 
"""[args.corrFullG, args.corrNormalG, args.corrWoTissuesG, args.corrTopG,
	args.corrTop100G, args.corrNRcvG]"""

### Get samples' correlations
s_corr = [args.corrFullS, args.corrNRcvS]
"""[args.corrFullS, args.corrNormalS, args.corrWoTissuesS, args.corrTopS,
	args.corrTop100S, args.corrNRcvS]"""

### Get pca's files
files_pca = [args.pcaFull, args.pcaNRcv]
"""[args.pcaRaw, args.pcaFiltered, args.pcaFull, args.pcaNormal,
	args.pcaWoTissues, args.pcaTop, args.pcaTop100, args.pcaNRcv]"""

### Get pca's images
images_pca = [args.IpcaFull, args.IpcaNRcv] 
"""[args.IpcaRaw, args.IpcaFiltered, args.IpcaFull, args.IpcaNormal,
	args.IpcaWoTissues, args.IpcaTop, args.IpcaTop100, args.IpcaNRcv]"""

### Get pca's HTML files
htmls_pca = [args.PpcaFull, args.PpcaNRcv]
"""[args.PpcaRaw, args.PpcaFiltered, args.PpcaFull, args.PpcaNormal,
	args.PpcaWoTissues, args.PpcaTop, args.PpcaTop100, args.PpcaNRcv]"""

### Get t-sne's files
files_tsne = [args.tsneFull, args.tsneNRcv]
"""[args.tsneRaw, args.tsneFiltered, args.tsneFull, args.tsneNormal,
	args.tsneWoTissues, args.tsneTop, args.tsneTop100, args.tsneNRcv]"""

### Get t-sne's images
images_tsne = [args.ItsneFull, args.ItsneNRcv]
"""[args.ItsneRaw, args.ItsneFiltered, args.ItsneFull, args.ItsneNormal, 
	args.ItsneWoTissues, args.ItsneTop, args.ItsneTop100, args.ItsneNRcv]"""

### Get t-sne's HTML files
htmls_tsne = [args.PtsneFull, args.PtsneNRcv]
"""[args.PtsneRaw, args.PtsneFiltered, args.PtsneFull, args.PtsneNormal, 
	args.PtsneWoTissues, args.PtsneTop, args.PtsneTop100, args.PtsneNRcv]"""

### Get UMAP's files
files_umap = [args.umapFull, args.umapNRcv]
"""[args.umapRaw, args.umapFiltered, args.umapFull, args.umapNormal,
	args.umapWoTissues, args.umapTop, args.umapTop100, args.umapNRcv]"""

### Get UMAP's images
images_umap = [args.IumapFull, args.IumapNRcv]
"""[args.IumapRaw, args.IumapFiltered, args.IumapFull, args.IumapNormal, 
	args.IumapWoTissues, args.IumapTop, args.IumapTop100, args.IumapNRcv]"""

### Get UMAP's HTML files
htmls_umap = [args.PumapFull, args.PumapNRcv]
"""[args.PumapRaw, args.PumapFiltered, args.PumapFull, args.PumapNormal, 
	args.PumapWoTissues, args.PumapTop, args.PumapTop100, args.PumapNRcv]"""

### Get QC images
#### General
general_qc_imgs = [args.IgeneralFull, args.IgeneralNRcv]
"""[args.IgenralRaw, args.IgenralFiltered, args.IgeneralFull, 
	args.IgeneralNormal, args.IgeneralWoTissues, args.IgeneralTop,
	args.IgeneralTop100, args.IgeneralNRcv]"""

#### Mean-variance
mv_imgs = [args.ImvFull, args.ImvNRcv]
"""[args.ImvRaw, args.ImvFiltered, args.ImvFull, 
	args.ImvNormal, args.ImvWoTissues, args.ImvTop,
	args.ImvTop100, args.ImvNRcv]"""

#### z_scores
zscores_imgs = [args.IzscoreFull, args.IzscoreNRcv]
"""[args.IzscoreRaw, args.IzscoreFiltered, args.IzscoreFull, 
	args.IzscoreNormal, args.IzscoreWoTissues, args.IzscoreTop, 
	args.IzscoreTop100, args.IzscoreNRcv]"""

### Get QC html files
#### General
general_qc_htmls = [args.PgeneralFull, args.PgeneralNRcv]
"""[args.PgenralRaw, args.PgenralFiltered, args.PgeneralFull, 
	args.PgeneralNormal, args.PgeneralWoTissues, args.PgeneralTop,
	args.PgeneralTop100, args.PgeneralNRcv]"""

#### Mean-variance
mv_htmls = [args.PmvFull, args.PmvNRcv]
"""[args.PmvRaw, args.PmvFiltered, args.PmvFull, 
	args.PmvNormal, args.PmvWoTissues, args.PmvTop,
	args.PmvTop100, args.PmvNRcv]"""

#### z_scores
zscores_htmls = [args.PzscoreFull, args.PzscoreNRcv]
"""[args.PzscoreRaw, args.PzscoreFiltered, args.PzscoreFull, 
	args.PzscoreNormal, args.PzscoreWoTissues, args.PzscoreTop, 
	args.PzscoreTop100, args.PzscoreNRcv]"""

### Get the clustermaps
s_clustermaps = [args.IclstrFullS, args.IclstrNRcvS]
"""[args.IclstrFullS, args.IclstrNormalS, args.IclstrWoTissueS, 
	args.IclstrTopS, args.IclstrTop100S, args.IclstrNRcvS]"""

g_clustermaps = [args.IclstrFullG, args.IclstrNRcvG]
"""[args.IclstrFullG, args.IclstrNormalG, args.IclstrWoTissueG, 
	args.IclstrTopG, args.IclstrTop100G, args.IclstrNRcvG]"""

sg_clustermaps = [args.IclstrFullSG, args.IclstrNRcvSG]
"""[args.IclstrFullSG, args.IclstrNormalSG, args.IclstrWoTissueSG, 
	args.IclstrTopSG, args.IclstrTop100SG, args.IclstrNRcvSG]"""