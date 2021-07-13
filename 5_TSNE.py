import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.manifold import TSNE 
from sklearn.preprocessing import StandardScaler

import Tools
import Config

def tsne(
    meta, rand, top100, cv_counts, raw_counts, top1000,
    counts_norm, filtered_counts, tissue_counts, p_tsne_cv,
    file_tsneCV, img_tsne_cv, p_tsne_raw, file_tsneRaw,
    img_tsne_raw, p_tsne_top, file_tsneTop, img_tsne_top,
    p_tsne_norm, file_tsneNorm, img_tsne_norm, normal, 
    p_tsne_normal, file_tsne_normal, img_tsne_normal, 
    p_tsne_rand, file_tsne_rand, img_tsne_rand, 
    p_tsne_top100, file_tsnetop100, img_tsne_top100, 
    p_tsne_tissues, file_tsne_tissues, img_tsne_tissues, 
    p_tsne_filtered, file_tsneFiltered, img_tsne_filtered, 
    counts_without_tissues, file_tsne_wo_tissues, 
    img_tsne_wo_tissues, p_tsne_wo_tissues):
    
    counts = [counts_norm, top1000, top100, cv_counts]
    #raw_counts, filtered_counts, normal, counts_without_tissues, 
    
    tissue_files = sorted([f for f in tissue_counts.iterdir() if f.is_file()])
    for path in tissue_files:
        counts.append(path)
    
    rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
    for path in rand_files:
        counts.append(path)

    tsvs = [file_tsneNorm, file_tsneTop, file_tsnetop100, file_tsneCV]
    #file_tsneRaw, file_tsneFiltered, file_tsne_normal, file_tsne_wo_tissues, 

    images = [img_tsne_norm, img_tsne_top, img_tsne_top100, img_tsne_cv]
    #img_tsne_raw, img_tsne_filtered, img_tsne_normal, img_tsne_wo_tissues, 

    htmls = [p_tsne_norm, p_tsne_top, p_tsne_top100, p_tsne_cv]
    #p_tsne_raw, p_tsne_filtered, p_tsne_normal, p_tsne_wo_tissues, 
    
    for i in range(len(tissue_files)):
        tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
        link_tsv = file_tsne_tissues.joinpath(tissue_name)
        link_img = img_tsne_tissues.joinpath(tissue_name)
        link_p = p_tsne_tissues.joinpath(tissue_name)
        
        Tools.create_folder(link_tsv)
        Tools.create_folder(link_img)
        Tools.create_folder(link_p)
        
        tsvs.append(link_tsv.joinpath(tissue_name + ".tsv"))
        images.append(link_img.joinpath(tissue_name + ".png"))
        htmls.append(link_p.joinpath(tissue_name + ".html"))

    for i in range(len(rand_files)):
        tsvs.append(file_tsne_rand.joinpath("random" + str(i) + ".tsv"))
        images.append(img_tsne_rand.joinpath("random" + str(i) + ".png"))
        htmls.append(p_tsne_rand.joinpath("random" + str(i) + ".html"))

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = '\t')
    tissues = metadata["smtsd"]

    for c, t, i, h in zip(counts, tsvs, images, htmls):
        f = pd.read_csv(c, header = 0, index_col = 0, sep = "\t")
        lib_size = f.sum(axis = 0).to_frame(name = "lib_size")
        f = pd.DataFrame(np.log2(f + 1))

        scaler = StandardScaler()
        std_counts = scaler.fit_transform(f.dropna().T)

        tsne = TSNE(n_components = 2)
        T = tsne.fit_transform(std_counts)
        df = pd.DataFrame(index = f.T.index)
        df["T1"] = T[:, 0]
        df["T2"] = T[:, 1]
        df = df.join(lib_size)
        df = df.join(tissues)
        df.sort_values("smtsd", inplace = True)

        fig = px.scatter(
            data_frame = df.dropna(), 
            x = "T1", y = "T2",
            color = df["smtsd"], 
            hover_data = [df.dropna().index, "lib_size"],
            #size = "lib_size",
            title = "GTEx t-sne")

        fig.write_html(str(h))
        fig.write_image(str(i), width = 2048, height = 1024)
        #fig.show()

        df.to_csv(t, sep="\t")

if __name__ == '__main__':	
    
    tsne(
        meta = Config.args.meta, rand = Config.args.rand,
        top100 = Config.args.top100, cv_counts = Config.args.cv,
        raw_counts = Config.args.bf, top1000 = Config.args.top1000,
        counts_norm = Config.args.norm,	filtered_counts = Config.args.af,
        tissue_counts = Config.args.tissue, p_tsne_cv = Config.args.PtsneCV,
        file_tsneCV = Config.args.tsneCV, img_tsne_cv = Config.args.ItsneCV,
        p_tsne_raw = Config.args.PtsneRaw, file_tsneRaw = Config.args.tsneRaw,
        img_tsne_raw = Config.args.ItsneRaw, p_tsne_top = Config.args.PtsneTop,
        file_tsneTop = Config.args.tsneTop, img_tsne_top = Config.args.ItsneTop,
        p_tsne_norm = Config.args.PtsneNorm, file_tsneNorm = Config.args.tsneNorm,
        img_tsne_norm = Config.args.ItsneNorm, normal = Config.args.onlyNormal, 
        p_tsne_normal = Config.args.PtsneNormal, file_tsne_normal = Config.args.tsneNormal, 
        img_tsne_normal = Config.args.ItsneNormal, p_tsne_rand = Config.args.PtsneRand,
        file_tsne_rand = Config.args.tsneRand, img_tsne_rand = Config.args.ItsneRand,
        p_tsne_top100 = Config.args.Ptsnetop100, file_tsnetop100 = Config.args.tsnetop100,
        img_tsne_top100 = Config.args.Itsnetop100, p_tsne_tissues = Config.args.PtsneTissue,
        file_tsne_tissues = Config.args.tsneTissue, img_tsne_tissues = Config.args.ItsneTissue,
        p_tsne_filtered = Config.args.PtsneFiltered, file_tsneFiltered = Config.args.tsneFiltered,
        img_tsne_filtered = Config.args.ItsneFiltered, counts_without_tissues = Config.args.WoTissues,
        file_tsne_wo_tissues = Config.args.tsneWoTissues, img_tsne_wo_tissues = Config.args.ItsneWoTissues,
        p_tsne_wo_tissues = Config.args.PtsneWoTissues)