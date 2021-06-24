import sys
import pandas as pd
import seaborn as  sns
import matplotlib.pyplot as plt

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
            bbox_to_anchor = (0, 1), loc = 'best', 
            bbox_transform = plt.gcf().transFigure)

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

        if c == 2:
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
            bbox_to_anchor = (0, 1), loc = 'best', 
            bbox_transform = plt.gcf().transFigure)

        g.savefig(str(i), dpi = 300)
        c += 1

def clustering_samples_genes(
    meta, cv_list, 
    counts_norm, top1000, cv_counts, rand, 
    sg_img_clstrRand, sg_img_clstrNorm, 
    sg_img_clstrTop, sg_img_clstrCv):
    
    counts = [counts_norm, top1000, cv_counts] 
    rand_files = sorted([f for f in g_corr_rand.iterdir() if f.is_file()])
    for path in rand_files:
        corr_matices.append(path)

    images = [sg_img_clstrNorm, sg_img_clstrTop, sg_img_clstrCv]
    for i in range(len(rand_files)):
        images.append(g_img_clstrRand.joinpath("random" + str(i) + ".png"))  

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = "\t")
    cv_list = pd.read_csv(cv_list, header = 0, index_col = 0, sep = "\t") 
    
    c = 0
    for m, i in zip(counts, images):
        f = pd.read_csv(m, header = 0, index_col = 0, sep = '\t')
        f = pd.DataFrame(np.log2(f + 1))

        count = f.T.join(metadata["smts"])
        count = count.dropna()
        tissues = count.pop("smts")

        palette1 = sns.hls_palette(10)
        palette2 = sns.color_palette("bwr",10)
        palette3 = sns.color_palette("inferno",10)
        palette = palette1 + palette2 + palette3
        lut = dict(zip(set(tissues.unique()), palette))
        col_colors = tissues.map(lut)

        if c == 2:
            count = count.T
            count = count.join(cv_list["GeneName"])
            count = count.join(cv_list["Class"])
            classes = count.pop("Class")

            paletteX = sns.color_palette("Set2" ,len(classes))
            lutX = dict(zip(set(classes.unique()), paletteX))
            row_colors = classes.map(lutX)

            yticklabels = count["GeneName"]
            data = count.iloc[:, count.columns != "GeneName"]

        else:
            row_colors, yticklabels = False, False
            data = count

        g = sns.clustermap(data, 
            vmin = max(count.max(axis = 1)), 
            vmax = min(count.min(axis = 1)), 
            row_colors = row_colors,
            col_colors = col_colors,
            cmap = "icefire",
            xticklabels = False, 
            yticklabels = yticklabels,
            method = "weighted",
            figsize = [15, 15])

        handles = [Patch(facecolor = lut[name]) for name in lut]
        g.ax_col_dendrogram.legend(handles, lut, title = 'Tissues',
            bbox_to_anchor = (0, 1), loc = 'best', 
            bbox_transform = plt.gcf().transFigure)

        if c == 2:
            handlesX = [Patch(facecolor = lutX[name]) for name in lutX]
            g.ax_row_dendrogram.legend(handlesX, lutX, title = 'Class',
                bbox_to_anchor = (0, 0), loc = 'best', 
                bbox_transform = plt.gcf().transFigure)

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

    clustering_samples_genes(
        counts_norm = Config.args.norm,
        top1000 = Config.args.top1000,
        cv_counts = Config.args.cv,
        rand = Config.args.rand,
        sg_img_clstrRand = Config.args.IclusterRandSG,
        sg_img_clstrNorm = Config.args.IclusterNormSG,
        sg_img_clstrTop = Config.args.IclusterTopSG,
        sg_img_clstrCv = Config.args.IclusterCVsg)