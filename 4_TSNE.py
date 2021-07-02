import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.manifold import TSNE 
from sklearn.preprocessing import StandardScaler

import Config

def tsne(
    meta, raw_counts, filtered_counts, counts_norm, top1000, cv_counts, 
    rand, file_tsneRaw, file_tsneFiltered, file_tsne_rand, file_tsneNorm, 
    file_tsneTop, file_tsneCV, img_tsne_raw, img_tsne_filtered, img_tsne_rand, 
    img_tsne_norm, img_tsne_top, img_tsne_cv, p_tsne_raw, p_tsne_filtered, 
    p_tsne_rand, p_tsne_norm, p_tsne_top, p_tsne_cv):
    
    counts = [raw_counts, filtered_counts, counts_norm, top1000, cv_counts]
    rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
    for path in rand_files:
        counts.append(path)

    tsvs = [file_tsneRaw, file_tsneFiltered, file_tsneNorm, file_tsneTop, file_tsneCV]
    images = [img_tsne_raw, img_tsne_filtered, img_tsne_norm, img_tsne_top, img_tsne_cv]
    htmls = [p_tsne_raw, p_tsne_filtered, p_tsne_norm, p_tsne_top, p_tsne_cv]

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
            color = "smtsd", 
            hover_data = [df.dropna().index, "lib_size"],
            #size = "lib_size",
            title = "GTEx t-sne")

        fig.write_html(str(h))
        fig.write_image(str(i), width = 2048, height = 1024)
        #fig.show()

        df.to_csv(t, sep="\t")

if __name__ == '__main__':	
    
    tsne(
        meta = Config.args.meta,
		raw_counts = Config.args.bf,
		filtered_counts = Config.args.af,
		counts_norm = Config.args.norm,
		top1000 = Config.args.top1000,
		cv_counts = Config.args.cv,
		rand = Config.args.rand,
        file_tsneRaw = Config.args.tsneRaw,
        file_tsneFiltered = Config.args.tsneFiltered,
        file_tsne_rand = Config.args.tsneRand,
        file_tsneNorm = Config.args.tsneNorm,
        file_tsneTop = Config.args.tsneTop,
        file_tsneCV = Config.args.tsneCV,
        img_tsne_raw = Config.args.ItsneRaw,
        img_tsne_filtered = Config.args.ItsneFiltered,
        img_tsne_rand = Config.args.ItsneRand,
        img_tsne_norm = Config.args.ItsneNorm,
        img_tsne_top = Config.args.ItsneTop,
        img_tsne_cv = Config.args.ItsneCV,
        p_tsne_raw = Config.args.PtsneRaw,
        p_tsne_filtered = Config.args.PtsneFiltered,
        p_tsne_rand = Config.args.PtsneRand,
        p_tsne_norm = Config.args.PtsneNorm,
        p_tsne_top = Config.args.PtsneTop,
        p_tsne_cv = Config.args.PtsneCV)