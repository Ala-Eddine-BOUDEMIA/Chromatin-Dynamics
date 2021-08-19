import gc
import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.manifold import TSNE 
from sklearn.preprocessing import StandardScaler

import Tools
import Config

def tsne(
    meta, counts, rand, tissues,
    p_tsne, img_tsne, files_tsne,
    p_tsne_rand, file_tsne_rand, 
    img_tsne_rand, p_tsne_tissues,
    file_tsne_tissues, img_tsne_tissues):

    for i, j, k in zip(img_tsne, p_tsne, files_tsne):
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
        link_tsv = file_tsne_tissues.joinpath(tissue_name)
        link_img = img_tsne_tissues.joinpath(tissue_name)
        link_p = p_tsne_tissues.joinpath(tissue_name)
        
        Tools.create_folder(link_tsv)
        Tools.create_folder(link_img)
        Tools.create_folder(link_p)
        
        files_tsne.append(link_tsv.joinpath(tissue_name + ".tsv"))
        img_tsne.append(link_img.joinpath(tissue_name + ".png"))
        p_tsne.append(link_p.joinpath(tissue_name + ".html"))

    Tools.create_folder(img_tsne_rand)
    Tools.create_folder(p_tsne_rand)
    Tools.create_folder(file_tsne_rand)
    for i in range(len(rand_files)):
        files_tsne.append(file_tsne_rand.joinpath("random" + str(i) + ".tsv"))
        img_tsne.append(img_tsne_rand.joinpath("random" + str(i) + ".png"))
        p_tsne.append(p_tsne_rand.joinpath("random" + str(i) + ".html"))

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = '\t')
    tissues = metadata[Config.args.smtsd]

    for c, t, i, h in zip(counts, files_tsne, img_tsne, p_tsne):
        f = pd.read_csv(c, header = 0, index_col = 0, sep = "\t")
        lib_size = f.sum(axis = 0).to_frame(name = "lib_size")
        f = pd.DataFrame(np.log2(f + 1))

        scaler = StandardScaler()
        std_counts = scaler.fit_transform(f.dropna().T)

        tsne = TSNE(n_components = 2, random_state = 0)
        T = tsne.fit_transform(std_counts)
        df = pd.DataFrame(index = f.T.index)
        df["T1"] = T[:, 0]
        df["T2"] = T[:, 1]
        df = df.join(lib_size)
        df = df.join(tissues)
        df.sort_values(Config.args.smtsd, inplace = True)

        fig = px.scatter(
            data_frame = df.dropna(), 
            x = "T1", y = "T2",
            color = df[Config.args.smtsd], 
            hover_data = [df.dropna().index, "lib_size"],
            #size = "lib_size",
            title = Config.args.dataset + " t-SNE")

        fig.write_html(str(h))
        fig.write_image(str(i), width = 2048, height = 1024)
        #fig.show()

        df.to_csv(t, sep="\t", float_format='%.3f')

        del(f)
        del(df)
        del(fig)
        del(std_counts)
        gc.collect()

if __name__ == '__main__':	
    
    tsne(
        meta = Config.args.meta, 
        counts = Config.counts,
        rand = Config.args.rand, 
        tissues = Config.args.tissue,
        p_tsne = Config.htmls_tsne,
        img_tsne = Config.images_tsne,
        files_tsne = Config.files_tsne,
        p_tsne_rand = Config.args.PtsneRand,
        file_tsne_rand = Config.args.tsneRand, 
        img_tsne_rand = Config.args.ItsneRand,
        p_tsne_tissues = Config.args.PtsneTissue,
        file_tsne_tissues = Config.args.tsneTissue, 
        img_tsne_tissues = Config.args.ItsneTissue)