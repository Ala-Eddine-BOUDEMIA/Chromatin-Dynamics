import gc
import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import Tools
import Config

def pca(
    meta, counts, rand, tissues,
    files_pca, images_pca, htmls_pca,
    p_pca_rand, img_pca_rand, file_pca_rand, 
    p_pca_tissue, img_pca_tissues, file_pca_tissue):

    for i, j, k in zip(images_pca, htmls_pca, files_pca):
        Tools.create_folder(('/').join(str(i)).split("/")[:-1])
        Tools.create_folder(('/').join(str(j)).split("/")[:-1])
        Tools.create_folder(('/').join(str(k)).split("/")[:-1])
    
    tissue_files = sorted([f for f in tissues.glob("**/*.tsv") if f.is_file()])
    for path in tissue_files:
        counts.append(path)
    
    rand_files = sorted([f for f in rand.glob("**/*.tsv") if f.is_file()])
    for path in rand_files:
        counts.append(path)

    for i in range(len(tissue_files)):
        tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
        link_tsv = file_pca_tissue.joinpath(tissue_name)
        link_img = img_pca_tissues.joinpath(tissue_name)
        link_p = p_pca_tissue.joinpath(tissue_name)
        
        Tools.create_folder(link_tsv)
        Tools.create_folder(link_img)
        Tools.create_folder(link_p)
        
        files_pca.append(link_tsv.joinpath("pca.tsv"))
        images_pca.append(link_img.joinpath("pca.png"))
        htmls_pca.append(link_p.joinpath("pca.html"))

    for i in range(len(rand_files)):
        link_tsv = file_pca_rand.joinpath("random" + str(i))
        link_img = img_pca_rand.joinpath("random" + str(i))
        link_html = p_pca_rand.joinpath("random" + str(i))

        Tools.create_folder(link_tsv)
        Tools.create_folder(link_img)
        Tools.create_folder(link_html)

        files_pca.append(link_tsv.joinpath("pca.tsv"))
        images_pca.append(link_img.joinpath("pca.png"))
        htmls_pca.append(link_html.joinpath("pca.html"))

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = '\t')
    sub_tissues = metadata[Config.args.smtsd]

    for c, t, i, h in zip(counts, files_pca, images_pca, htmls_pca):
        f = pd.read_csv(c, header = 0, index_col = 0, sep = "\t")
        lib_size = f.sum(axis = 0).to_frame(name = "lib_size")
        f = pd.DataFrame(np.log2(f + 1))
        scaler = StandardScaler()
        std_counts = scaler.fit_transform(f.T.dropna())

        pca = PCA(n_components = 2)
        P = pca.fit_transform(std_counts)
        ratio = pca.explained_variance_ratio_ * 100

        d = pd.DataFrame(index = f.T.dropna().index)
        d["PC1"] = P[:, 0]
        d["PC2"] = P[:, 1]
        d = d.join(lib_size)
        d = d.join(sub_tissues)
        d.sort_values(Config.args.smtsd, inplace = True)

        print("PC1: ", round(ratio[0], 2))
        print("PC2: ", round(ratio[1], 2))

        loading_scores = pd.Series(pca.components_[0], index = f.index)
        sorted_loading_scores = loading_scores.abs().sort_values(ascending = False)
        top_10_genes = sorted_loading_scores[0:10].index.values
        print(loading_scores[top_10_genes])
        
        fig = px.scatter(
            data_frame = d.dropna(), 
            x = "PC1", y = "PC2",
            color = d[Config.args.smtsd], 
            hover_data = [d.dropna().index, "lib_size"],
            #size = "lib_size",
            title = Config.args.dataset + "  PCA")

        fig.write_html(str(h))
        fig.write_image(str(i), width = 2048, height = 1024)
        #fig.show()

        d.to_csv(str(t), sep="\t", float_format='%.3f')

        del(f)
        del(d)
        del(fig)
        del(std_counts)
        gc.collect()

if __name__ == '__main__':	
    
    pca(
        meta = Config.args.meta,
        counts = Config.counts, 
        rand = Config.args.rand,
        tissues = Config.args.tissue,
        files_pca = Config.files_pca,
        images_pca = Config.images_pca, 
        htmls_pca = Config.htmls_pca,
        p_pca_rand = Config.args.PpcaRand,
        img_pca_rand = Config.args.IpcaRand, 
        file_pca_rand = Config.args.pcaRand, 
        p_pca_tissue = Config.args.PpcaTissue,
        img_pca_tissues = Config.args.IpcaTissue,
        file_pca_tissue = Config.args.pcaTissue) 