import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import Tools
import Config

def pca(
    meta, raw_counts, filtered_counts, counts_norm, top1000, cv_counts, rand, file_pcaRaw,
    file_pcaFiltered, file_pca_rand, file_pcaNorm, file_pcaTop, file_pcaCV, img_pca_raw, img_pca_filtered, 
    img_pca_rand, img_pca_norm, img_pca_top, img_pca_cv, p_pca_raw, p_pca_filtered, p_pca_rand, p_pca_norm,
    p_pca_top, p_pca_cv):
    
    counts = [raw_counts, filtered_counts, counts_norm, top1000, cv_counts]
    rand_files = Tools.parse_dir(rand)
    counts.append(rand_files)

    tsvs = [file_pcaRaw, file_pcaFiltered, file_pcaNorm, file_pcaTop, file_pcaCV]
    images = [img_pca_raw, img_pca_filtered, img_pca_norm, img_pca_top, img_pca_cv]
    htmls = [p_pca_raw, p_pca_filtered, p_pca_norm, p_pca_top, p_pca_cv]

    for i in range(len(rand_files)):
        tsvs.append(file_pca_rand.joinpath("random" + str(i) + ".tsv"))
        images.append(img_pca_rand.joinpath("random" + str(i) + ".png"))
        htmls.append(p_pca_rand.joinpath("random" + str(i)) + ".html")

    metadata = pd.read_csv(meta, sep = '\t')
    tissues = metadata["smtsd"]

    for c, t, i, h in zip(counts, tsvs, images, htmls):
        f = pd.read_csv(c, sep = "\t")
        
        lib_size = f.sum(axis = 0).to_frame(name = "lib_size")
        scaler = StandardScaler()
        std_counts = scaler.fit_transform(f.dropna().T)

        pca = PCA(n_components = 2)
        P = pca.fit_transform(std_counts)
        ratio = pca.explained_variance_ratio_ * 100

        d = pd.DataFrame(index = f.T.index)
        d["PC1"] = P[:, 0]
        d["PC2"] = P[:, 1]
        d = d.join(lib_size)
        d = d.join(tissues)
        d.sort_values("smtsd", inplace=True)

        print("PC1: ", round(ratio[0],2))
        print("PC2: ", round(ratio[1],2))

        loading_scores = pd.Series(pca.components_[0], index = f.index)
        sorted_loading_scores = loading_scores.abs().sort_values(ascending = False)
        top_10_genes = sorted_loading_scores[0:10].index.values
        print(loading_scores[top_10_genes])
        fig = px.scatter(
            d.dropna(), 
            x="PC1", y="PC2",
            color="smtsd", 
            hover_data=[d.dropna().index, "lib_size"],
            #size="lib_size",
            title="GTEx PCA")

        fig.write_html(str(h))
        fig.write_image(str(i))
        fig.show()

        d.to_csv(t, sep="\t")

if __name__ == '__main__':	
    
    pca(
        meta = Config.args.meta,
		raw_counts = Config.args.bf,
		filtered_counts = Config.args.af,
		counts_norm = Config.args.norm,
		top1000 = Config.args.top1000,
		cv_counts = Config.args.cv,
		rand = Config.args.rand,
        file_pcaRaw = Config.args.pcaRaw,
        file_pcaFiltered = Config.args.pcaFiltered,
        file_pca_rand = Config.args.pcaRand,
        file_pcaNorm = Config.args.pcaNorm,
        file_pcaTop = Config.args.pcaTop,
        file_pcaCV = Config.args.pcaCV,
        img_pca_raw = Config.args.IpcaRaw,
        img_pca_filtered = Config.args.IpcaFiltered,
        img_pca_rand = Config.args.IpcaRand,
        img_pca_norm = Config.args.IpcaNorm,
        img_pca_top = Config.args.IpcaTop,
        img_pca_cv = Config.args.IpcaCV,
        p_pca_raw = Config.args.PpcaRaw,
        p_pca_filtered = Config.args.PpcaFiltered,
        p_pca_rand = Config.args.PpcaRand,
        p_pca_norm = Config.args.PpcaNorm,
        p_pca_top = Config.args.PpcaTop,
        p_pca_cv = Config.args.PpcaCV)