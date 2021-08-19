import gc
import numpy as np
import pandas as pd
import plotly.express as px

from umap import UMAP
from sklearn.preprocessing import StandardScaler

import Tools
import Config

def umap(
    meta, counts, rand, tissues,
    p_umap, img_umap, files_umap,
    p_umap_rand, file_umap_rand, 
    img_umap_rand, p_umap_tissues,
    file_umap_tissues, img_umap_tissues):

    for i, j, k in zip(img_umap, p_umap, files_umap):
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
        link_tsv = file_umap_tissues.joinpath(tissue_name)
        link_img = img_umap_tissues.joinpath(tissue_name)
        link_p = p_umap_tissues.joinpath(tissue_name)
        
        Tools.create_folder(link_tsv)
        Tools.create_folder(link_img)
        Tools.create_folder(link_p)
        
        files_umap.append(link_tsv.joinpath(tissue_name + ".tsv"))
        img_umap.append(link_img.joinpath(tissue_name + ".png"))
        p_umap.append(link_p.joinpath(tissue_name + ".html"))

    Tools.create_folder(img_umap_rand)
    Tools.create_folder(p_umap_rand)
    Tools.create_folder(file_umap_rand)
    for i in range(len(rand_files)):
        files_umap.append(file_umap_rand.joinpath("random" + str(i) + ".tsv"))
        img_umap.append(img_umap_rand.joinpath("random" + str(i) + ".png"))
        p_umap.append(p_umap_rand.joinpath("random" + str(i) + ".html"))

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = '\t')
    tissues = metadata[Config.args.smtsd]

    for c, t, i, h in zip(counts, files_umap, img_umap, p_umap):
        print(c)
        f = pd.read_csv(c, header = 0, index_col = 0, sep = "\t")
        lib_size = f.sum(axis = 0).to_frame(name = "lib_size")
        f = pd.DataFrame(np.log2(f + 1))

        scaler = StandardScaler()
        std_counts = scaler.fit_transform(f.dropna().T)
        
        reducer = UMAP(n_components = 2, random_state = 0)
        U = reducer.fit_transform(std_counts)
        df = pd.DataFrame(index = f.T.index)
        df["U1"] = U[:, 0]
        df["U2"] = U[:, 1]
        df = df.join(lib_size)
        df = df.join(tissues)
        df = df.join(metadata['database'])
        df.sort_values(Config.args.smtsd, inplace = True)

        fig = px.scatter(
            data_frame = df.dropna(), 
            x = "U1", y = "U2",
            color = "database", 
            hover_data = [df.dropna().index, "lib_size"],
            #size = "lib_size",
            title = Config.args.dataset + " umap")

        #fig.write_html(str(h))
        #fig.write_image(str(i), width = 2048, height = 1024)
        fig.show()

        #df.to_csv(t, sep="\t", float_format='%.3f')

        del(f)
        del(df)
        del(fig)
        del(std_counts)
        gc.collect()

if __name__ == '__main__':	
    
    umap(
        meta = Config.args.meta, 
        counts = Config.counts,
        rand = Config.args.rand, 
        tissues = Config.args.tissue,
        p_umap = Config.htmls_umap,
        img_umap = Config.images_umap,
        files_umap = Config.files_umap,
        p_umap_rand = Config.args.PumapRand,
        file_umap_rand = Config.args.umapRand, 
        img_umap_rand = Config.args.IumapRand,
        p_umap_tissues = Config.args.PumapTissue,
        file_umap_tissues = Config.args.umapTissue, 
        img_umap_tissues = Config.args.IumapTissue)