import sys
import pandas as pd
import seaborn as  sns

from matplotlib.patches import Patch

import Config

sys.setrecursionlimit(1000000)

def clustering_samples(
    meta, s_corr_norm, s_corr_top1000, s_corr_rand, s_corr_cv,
    s_img_clstrRand, s_img_clstrNorm, s_img_clstrTop, s_img_clstrCv):
    
    corr_matices = [s_corr_norm, s_corr_top1000, s_corr_cv] 
    rand_files = sorted([f for f in s_corr_rand.iterdir() if f.is_file()])
    for path in rand_files:
        corr_matices.append(path)

    images = [s_img_clstrNorm, s_img_clstrTop, s_img_clstrCv]  
    for i in range(len(rand_files)):
        images.append(s_img_clstrRand.joinpath("random" + str(i) + ".png"))

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = "\t")
    
    for m, i in zip(corr_matices, images):
        f = pd.read_csv(m, header = 0, index_col = 0, sep = '\t')

        correlation_matrix = f.join(metadata["smts"])
        correlation_matrix = correlation_matrix.dropna()
        tissues = correlation_matrix.pop("smts")

        palette1 = sns.hls_palette(10)
        palette2 = sns.color_palette("bwr",10)
        palette3 = sns.color_palette("inferno",10)
        palette = palette1 + palette2 + palette3
        lut = dict(zip(set(tissues.unique()), palette))
        colors = tissues.map(lut)

        g = sns.clustermap(correlation_matrix, 
            vmin = max(correlation_matrix.max(axis = 1)), 
            vmax = min(correlation_matrix.min(axis = 1)), 
            cmap = "icefire", standard_scale = 1,
            row_colors = colors, col_colors = colors, 
            xticklabels = False, yticklabels = False,
            method = "single")
        
        handles = [Patch(facecolor = lut[name]) for name in lut]
        g.ax_row_dendrogram.legend(handles, lut, title = 'Tissues',
            bbox_to_anchor = (0, 1), loc='best')

        g.savefig(str(i), dpi = 300)

def clustering_genes(
    cv_list, g_corr_norm, g_corr_top1000, g_corr_cv, g_corr_rand,
    g_img_clstrRand, g_img_clstrNorm, g_img_clstrTop, g_img_clstrCv):
    
    corr_matices = [g_corr_norm, g_corr_top1000, g_corr_cv] 
    rand_files = sorted([f for f in g_corr_rand.iterdir() if f.is_file()])
    for path in rand_files:
        corr_matices.append(path)

    images = [g_img_clstrNorm, g_img_clstrTop, g_img_clstrCv]
    for i in range(len(rand_files)):
        images.append(g_img_clstrRand.joinpath("random" + str(i) + ".png"))

    metadata = pd.read_csv(cv_list, header = 0, index_col = 0, sep = "\t")
    
    c = 0
    for m, i in zip(corr_matices, images):
        f = pd.read_csv(m, header = 0, index_col = 0, sep = '\t')

        if c < 3:
            correlation_matrix = f.join(metadata["Class"])
            correlation_matrix = f.join(metadata["GeneName"])
            correlation_matrix = correlation_matrix.dropna()
            classes = correlation_matrix.pop("Class")

            palette = sns.hls_palette(len(classes))
            lut = dict(zip(set(classes.unique()), palette))
            colors = classes.map(lut)
            labels = correlation_matrix["GeneName"]

            data = correlation_matrix.iloc[:, correlation_matrix.columns != "GeneName"]
        
        else :
            colors, labels = False,  False
            data = correlation_matrix

        g = sns.clustermap(data, 
            vmin = max(correlation_matrix.max(axis = 1)), 
            vmax = min(correlation_matrix.min(axis = 1)), 
            cmap = "icefire", standard_scale = 1,
            row_colors = colors, col_colors = colors, 
            xticklabels = labels, yticklabels = labels,
            method = "single")
        
        handles = [Patch(facecolor = lut[name]) for name in lut]
        g.ax_row_dendrogram.legend(handles, lut, title = 'Class',
            bbox_to_anchor = (0, 1), loc='best')

        g.savefig(str(i), dpi = 300)
        c += 1

if __name__ == '__main__':

    clustering_samples(
        meta = Config.args.meta,
        s_corr_norm = Config.args.corrNormS,
		s_corr_top1000 = Config.args.corrTopS,
		s_corr_rand = Config.args.corrRandS,
		s_corr_cv = Config.args.corrCVs,
        s_img_clstrRand = Config.args.IclusterRandS,
        s_img_clstrNorm = Config.args.IclusterNormS,
        s_img_clstrTop = Config.args.IclusterTopS,
        s_img_clstrCv = Config.args.IclusterCVs)

    clustering_genes(
        cv_list = Config.args.list,
        g_corr_norm = Config.args.corrNormG,
		g_corr_top1000 = Config.args.corrTopG,
		g_corr_cv = Config.args.corrCVg,
		g_corr_rand = Config.args.corrRandG,
        g_img_clstrRand = Config.args.IclusterRandG,
        g_img_clstrNorm = Config.args.IclusterNormG,
        g_img_clstrTop = Config.args.IclusterTopG,
        g_img_clstrCv = Config.args.IclusterCVG)