import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import Tools
import Config

def pca(
    meta, rand, top88, cv_counts, raw_counts, top1000, 
    counts_norm, filtered_counts, tissue_counts, p_pca_cv,
    file_pcaCV, img_pca_cv, p_pca_raw, file_pcaRaw, img_pca_raw, 
    p_pca_top, file_pcaTop, img_pca_top, p_pca_norm, file_pcaNorm, 
    img_pca_norm, p_pca_rand, file_pca_rand, img_pca_rand, p_pca_top88, 
    file_pcaTop88, img_pca_top88, file_pca_tissue, p_pca_tissue, 
    img_pca_tissues, p_pca_filtered, file_pcaFiltered, img_pca_filtered):
    
    counts = [raw_counts, filtered_counts, counts_norm, 
        top1000, top88, cv_counts]

    tissue_files = sorted([f for f in tissue_counts.iterdir() if f.is_file()])
    for path in tissue_files:
        counts.append(path)

    rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
    for path in rand_files:
        counts.append(path)

    tsvs = [file_pcaRaw, file_pcaFiltered, file_pcaNorm, 
        file_pcaTop, file_pcaTop88, file_pcaCV]
    
    images = [img_pca_raw, img_pca_filtered, img_pca_norm, 
        img_pca_top, img_pca_top88, img_pca_cv]
    
    htmls = [p_pca_raw, p_pca_filtered, p_pca_norm, 
        p_pca_top, p_pca_top88, p_pca_cv]

    for i in range(len(tissue_files)):
        tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
        link_tsv = file_pca_tissue.joinpath(tissue_name)
        link_img = img_pca_tissues.joinpath(tissue_name)
        link_p = p_pca_tissue.joinpath(tissue_name)
        
        Tools.create_folder(link_tsv)
        Tools.create_folder(link_img)
        Tools.create_folder(link_p)
        
        tsvs.append(link_tsv)
        images.append(link_img)
        htmls.append(link_p)

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
    sub_tissues = metadata["smtsd"]

    for c, t, i, h in zip(counts, tsvs, images, htmls):
        f = pd.read_csv(c, header = 0, index_col = 0, sep = "\t")
        lib_size = f.sum(axis = 0).to_frame(name = "lib_size")
        f = pd.DataFrame(np.log2(f + 1)).T

        scaler = StandardScaler()
        std_counts = scaler.fit_transform(f.T.dropna())

        pca = PCA(n_components = 2)
        P = pca.fit_transform(std_counts)
        ratio = pca.explained_variance_ratio_ * 100

        d = pd.DataFrame(index = f.T.index)
        d["PC1"] = P[:, 0]
        d["PC2"] = P[:, 1]
        d = d.join(lib_size)
        d = d.join(sub_tissues)
        d.sort_values("smtsd", inplace = True)

        print("PC1: ", round(ratio[0], 2))
        print("PC2: ", round(ratio[1], 2))

        loading_scores = pd.Series(pca.components_[0], index = f.index)
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

        fig.write_html(h.joinpath("pca.html"))
        fig.write_image(i.joinpath("pca.png"), width = 2048, height = 1024)
        fig.show()

        d.to_csv(t.joinpath("pca.tsv"), sep="\t")

if __name__ == '__main__':	
    
    pca(
        meta = Config.args.meta, rand = Config.args.rand,
        top88 = Config.args.top88, cv_counts = Config.args.cv,
        raw_counts = Config.args.bf, top1000 = Config.args.top1000, 
        counts_norm = Config.args.norm, filtered_counts = Config.args.af, 
        tissue_counts = Config.args.tissue, p_pca_cv = Config.args.PpcaCV,
        file_pcaCV = Config.args.pcaCV, img_pca_cv = Config.args.IpcaCV,
        p_pca_raw = Config.args.PpcaRaw, file_pcaRaw = Config.args.pcaRaw, 
        img_pca_raw = Config.args.IpcaRaw, p_pca_top = Config.args.PpcaTop,
        file_pcaTop = Config.args.pcaTop, img_pca_top = Config.args.IpcaTop,
        p_pca_norm = Config.args.PpcaNorm, file_pcaNorm = Config.args.pcaNorm, 
        img_pca_norm = Config.args.IpcaNorm, p_pca_rand = Config.args.PpcaRand,
        file_pca_rand = Config.args.pcaRand, img_pca_rand = Config.args.IpcaRand, 
        p_pca_top88 = Config.args.PpcaTop88, file_pcaTop88 = Config.args.pcaTop88, 
        img_pca_top88 = Config.args.IpcaTop88, p_pca_tissue = Config.args.PpcaTissue,
        file_pca_tissues = Config.args.pcaTissue, img_pca_tissues = Config.args.IpcaTissue,
        p_pca_filtered = Config.args.PpcaFiltered, file_pcaFiltered = Config.args.pcaFiltered,
        img_pca_filtered = Config.args.IpcaFiltered)