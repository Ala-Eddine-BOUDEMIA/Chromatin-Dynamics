import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import Tools
import Config

def pca(
    meta, raw_counts, filtered_counts, counts_norm, top1000, cv_counts, rand, 
    file_pcaRaw, file_pcaFiltered, file_pca_rand, file_pcaNorm, file_pcaTop, 
    file_pcaCV, img_pca_raw, img_pca_filtered, img_pca_rand, img_pca_norm, 
    img_pca_top, img_pca_cv, p_pca_raw, p_pca_filtered, p_pca_rand, p_pca_norm, 
    p_pca_top, p_pca_cv):
    
    counts = [raw_counts, filtered_counts, counts_norm, top1000, cv_counts]
    rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
    for path in rand_files:
        counts.append(path)

    tsvs = [file_pcaRaw, file_pcaFiltered, file_pcaNorm, file_pcaTop, file_pcaCV]
    images = [img_pca_raw, img_pca_filtered, img_pca_norm, img_pca_top, img_pca_cv]
    htmls = [p_pca_raw, p_pca_filtered, p_pca_norm, p_pca_top, p_pca_cv]

    for i in range(len(rand_files)):
        link_tsv = file_pca_rand.joinpath("random" + str(i) + "/")
        link_img = img_pca_rand.joinpath("random" + str(i) + "/")
        link_html = p_pca_rand.joinpath("random" + str(i) + "/")

        Tools.create_folder(link_tsv)
        Tools.create_folder(link_img)
        Tools.create_folder(link_html)

        tsvs.append(link_tsv)
        images.append(link_img)
        htmls.append(link_html)

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = '\t')
    tissues = metadata["smts"]
    sub_tissues = metadata["smtsd"]

    for c, tsv, i, h in zip(counts, tsvs, images, htmls):
        f = pd.read_csv(c, header = 0, index_col = 0, sep = "\t")
        lib_size = f.sum(axis = 0).to_frame(name = "lib_size")
        f = pd.DataFrame(np.log2(f + 1)).T
        f = f.join(tissues)

        tissue_types = tissues.unique()
        for t in tissue_types:
            df = pd.DataFrame(columns = f.columns)

            for i in f.index:
                if f.loc[i]["smts"] == t:
                    df = df.append(f.loc[i])

            df.pop("smts")

            scaler = StandardScaler()
            std_counts = scaler.fit_transform(df.dropna())

            pca = PCA(n_components = 2)
            P = pca.fit_transform(std_counts)
            ratio = pca.explained_variance_ratio_ * 100

            d = pd.DataFrame(index = df.index)
            d["PC1"] = P[:, 0]
            d["PC2"] = P[:, 1]
            d = d.join(lib_size)
            d = d.join(sub_tissues)
            d.sort_values("smtsd", inplace = True)

            print("PC1: ", round(ratio[0], 2))
            print("PC2: ", round(ratio[1], 2))

            loading_scores = pd.Series(pca.components_[0], index = df.index)
            sorted_loading_scores = loading_scores.abs().sort_values(ascending = False)
            top_10_genes = sorted_loading_scores[0:10].index.values
            print(loading_scores[top_10_genes])
            
            fig = px.scatter(
                data_frame = d.dropna(), 
                x = "PC1", y = "PC2",
                color = "smtsd", 
                hover_data = [d.dropna().index, "lib_size"],
                #size = "lib_size",
                title = "GTEx PCA")

            fig.write_html(h.joinpath("pca_" + t + ".html"))
            fig.write_image(i.joinpath("pca_" + t + ".png"), width = 2048, height = 1024)
            fig.show()

            d.to_csv(tsv.joinpath("pca_" + t + ".tsv"), sep="\t")

if __name__ == '__main__':	
    
    pca(
        meta = Config.args.meta, raw_counts = Config.args.bf,
		filtered_counts = Config.args.af, counts_norm = Config.args.norm,
		top1000 = Config.args.top1000, cv_counts = Config.args.cv, rand = Config.args.rand,
        file_pcaRaw = Config.args.pcaRaw, file_pcaFiltered = Config.args.pcaFiltered,
        file_pca_rand = Config.args.pcaRand, file_pcaNorm = Config.args.pcaNorm,
        file_pcaTop = Config.args.pcaTop, file_pcaCV = Config.args.pcaCV,
        img_pca_raw = Config.args.IpcaRaw, img_pca_filtered = Config.args.IpcaFiltered,
        img_pca_rand = Config.args.IpcaRand, img_pca_norm = Config.args.IpcaNorm,
        img_pca_top = Config.args.IpcaTop, img_pca_cv = Config.args.IpcaCV,
        p_pca_raw = Config.args.PpcaRaw, p_pca_filtered = Config.args.PpcaFiltered,
        p_pca_rand = Config.args.PpcaRand, p_pca_norm = Config.args.PpcaNorm,
        p_pca_top = Config.args.PpcaTop, p_pca_cv = Config.args.PpcaCV)