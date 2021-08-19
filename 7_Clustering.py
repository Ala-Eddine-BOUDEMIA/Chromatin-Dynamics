import gc
import sys
import numpy as np
import pandas as pd
import seaborn as  sns
import matplotlib.pyplot as plt

from matplotlib.patches import Patch
from sklearn.preprocessing import StandardScaler

import Tools
import Config

sys.setrecursionlimit(1000000)

def clustering_samples(
    meta, s_corr, s_corr_rand, s_corr_tissues,
    s_clustermaps, s_img_clstrRand, s_img_clstrTissues):

    for i in s_clustermaps:
        Tools.create_folder(('/').join(str(i)).split("/")[:-1])

    tissue_files = sorted([f for f in s_corr_tissues.glob('**/*.tsv') if f.is_file()])
    for path in tissue_files:
        s_corr.append(path) 
    
    rand_files = sorted([f for f in s_corr_rand.glob('**/*.tsv') if f.is_file()])
    for path in rand_files:
        s_corr.append(path)

    for i in range(len(tissue_files)):
        tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
        link = s_img_clstrTissues.joinpath(tissue_name)
        Tools.create_folder(link)
        s_clustermaps.append(link.joinpath(tissue_name + ".png"))

    Tools.create_folder(s_img_clstrRand)
    for i in range(len(rand_files)):
        s_clustermaps.append(s_img_clstrRand.joinpath("random" + str(i) + ".png"))

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = "\t")
    c = 0
    for m, i in zip(s_corr, s_clustermaps):
        correlation_matrix  = pd.read_csv(m, header = 0, index_col = 0, sep = '\t')

        if c <= 2: 
            tissue_type = Config.args.smts
        else:
            tissue_type = Config.args.smtsd

        correlation_matrix = correlation_matrix.join(metadata[tissue_type])
        correlation_matrix = correlation_matrix.dropna()
        tissues = correlation_matrix.pop(tissue_type)

        palette1 = sns.hls_palette(10)
        palette2 = sns.color_palette("bwr", 4)
        palette3 = sns.color_palette("Paired", 12)
        palette4 = sns.color_palette("Set2",10)

        palette = palette1 + palette2 + palette3 + palette4 
        lut = dict(zip(set(tissues.unique()), palette))
        colors = tissues.map(lut)

        g = sns.clustermap(correlation_matrix, 
            vmin = -1, 
            vmax = 1, 
            cmap = "icefire", metric = Config.distance_metric,
            row_colors = colors, col_colors = colors, 
            xticklabels = False, yticklabels = False,
            method = "average")
        
        handles = [Patch(facecolor = lut[name]) for name in lut]
        g.ax_row_dendrogram.legend(handles, lut, title = 'Tissues',
            bbox_to_anchor = (0, 1), loc = 'best', 
            bbox_transform = plt.gcf().transFigure)

        g.savefig(str(i), dpi = 300)
        plt.close('all')

        # Free memory
        del(correlation_matrix)
        del(g)
        gc.collect()

def clustering_genes(
    cv_list, g_corr, g_corr_rand, g_corr_tissue,
    g_clustermaps, g_img_clstrRand, g_img_clstrTissues):

    for i in g_clustermaps:
        Tools.create_folder(('/').join(str(i)).split("/")[:-1])

    tissue_files = sorted([f for f in g_corr_tissue.glob('**/*.tsv') if f.is_file()])
    for path in tissue_files:
        g_corr.append(path)
    
    rand_files = sorted([f for f in g_corr_rand.glob('**/*.tsv') if f.is_file()])
    for path in rand_files:
        g_corr.append(path)
    
    for i in range(len(tissue_files)):
        tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
        link = g_img_clstrTissues.joinpath(tissue_name)
        Tools.create_folder(link)
        g_clustermaps.append(link.joinpath(tissue_name + ".png"))
    
    Tools.create_folder(g_img_clstrRand)
    for i in range(len(rand_files)):
        g_clustermaps.append(g_img_clstrRand.joinpath("random" + str(i) + ".png"))

    metadata = pd.read_csv(cv_list, 
        header = 0, index_col = 0, sep = ";")

    c = 0
    for m, i in zip(g_corr, g_clustermaps):
        correlation_matrix = pd.read_csv(m, header = 0, index_col = 0, sep = '\t')
        # change c values depending on the config file
        if c == 1 or c == 2 or c > 36:
            g = sns.clustermap(correlation_matrix, 
                vmin = -1, 
                vmax = 1, 
                cmap = sns.color_palette("vlag", as_cmap = True), 
                metric = Config.distance_metric,
                xticklabels = False, yticklabels = False,
                method = "average", figsize = [25, 25])
        else:
            correlation_matrix = correlation_matrix.join(metadata["Class"])
            correlation_matrix = correlation_matrix.join(metadata["GeneName"])
            correlation_matrix = correlation_matrix.dropna()
            classes = correlation_matrix.pop("Class")

            palette = sns.color_palette("Set2" ,len(classes))
            lut = dict(zip(set(classes.unique()), palette))
            colors = classes.map(lut)
            labels = correlation_matrix["GeneName"]

            data = correlation_matrix.iloc[:, correlation_matrix.columns != "GeneName"]

            g = sns.clustermap(data, 
                vmin = -1, 
                vmax = 1, 
                cmap = sns.color_palette("vlag", as_cmap = True), 
                metric = Config.distance_metric,
                row_colors = colors, col_colors = colors, 
                xticklabels = labels, yticklabels = labels,
                method = "average", figsize = [25, 25])

            handles = [Patch(facecolor = lut[name]) for name in lut]
            g.ax_row_dendrogram.legend(handles, lut, title = 'Class',
                bbox_to_anchor = (0, 1), loc = 'best', 
                bbox_transform = plt.gcf().transFigure)

        g.savefig(str(i), dpi = 300)
        plt.close('all')
        c += 1

        # Free memory
        del(correlation_matrix)
        del(g)
        gc.collect()

def clustering_samples_genes(
    meta, cv_list, counts, rand, by_tissue,
    sg_clustermaps, sg_img_clstrRand, sg_img_clstrTissues):

    for i in sg_clustermaps:
        Tools.create_folder(('/').join(str(i)).split("/")[:-1])

    tissue_files = sorted([f for f in by_tissue.glob('**/*.tsv') if f.is_file()])
    for path in tissue_files:
        counts.append(path)
    
    rand_files = sorted([f for f in rand.glob('**/*.tsv') if f.is_file()])
    for path in rand_files:
        counts.append(path)
    
    for i in range(len(tissue_files)):
        tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
        link = sg_img_clstrTissues.joinpath(tissue_name)
        Tools.create_folder(link)
        sg_clustermaps.append(link.joinpath(tissue_name + ".png"))
    
    Tools.create_folder(sg_img_clstrRand)
    for i in range(len(rand_files)):
        sg_clustermaps.append(sg_img_clstrRand.joinpath("random" + str(i) + ".png"))

    metadata = pd.read_csv(meta, header = 0, index_col = 0, sep = "\t")
    cv_list = pd.read_csv(cv_list, header = 0, index_col = 0, sep = ";") 
    
    c = 0
    for m, i in zip(counts, sg_clustermaps):
        count = pd.read_csv(m, header = 0, index_col = 0, sep = '\t')
        count = pd.DataFrame(np.log2(count + 1))

        count = count.T

        if c <= 2 or c > 36:
            tissue_type = Config.args.smts
        else:
            tissue_type = Config.args.smtsd
        
        count = count.join(metadata[tissue_type])
        count = count.join(metadata["gdc_cases.samples.sample_type"])
        count = count.dropna()
        tissues = count.pop(tissue_type)
        samples = count.pop("gdc_cases.samples.sample_type")

        palette1 = sns.hls_palette(10)
        palette2 = sns.color_palette("bwr", 4)
        palette3 = sns.color_palette("Paired", 12)
        palette4 = sns.color_palette("Set2",10)
        palette = palette1 + palette2 + palette3 + palette4
        lut = dict(zip(set(tissues.unique()), palette))
        col_colors1 = tissues.map(lut)

        paletteY = sns.color_palette("Set2",10)
        lutY = dict(zip(set(samples.unique()), paletteY))
        col_colors2 = samples.map(lutY)

        col_colors = [col_colors1, col_colors2]
        
        # change c values depending on the config file
        if c > 36:
            g = sns.clustermap(count.T, 
                vmin = min(data.max(axis = 1)), 
                vmax = max(data.min(axis = 1)),  
                col_colors = col_colors,
                cmap = "icefire", 
                metric = Config.distance_metric,
                xticklabels = False, yticklabels = False,
                method = "average", figsize = [25, 25])
        else:
            count = count.T
            count = count.join(cv_list["GeneName"])
            count = count.join(cv_list["Class"])
            classes = count.pop("Class")

            paletteX = sns.color_palette("Set2" ,len(classes))
            lutX = dict(zip(set(classes.unique()), paletteX))
            row_colors = classes.map(lutX)

            yticklabels = count["GeneName"]
            data = count.iloc[:, count.columns != "GeneName"]

            g = sns.clustermap(data, 
                vmin = max(data.max(axis = 1)), 
                vmax = min(data.min(axis = 1)), 
                row_colors = row_colors,
                col_colors = col_colors,
                cmap = "icefire",
                metric = Config.distance_metric,
                xticklabels = False, 
                yticklabels = yticklabels,
                method = "average",
                figsize = [25, 25])

            handlesX = [Patch(facecolor = lutX[name]) for name in lutX]
            g.ax_row_dendrogram.legend(handlesX, lutX, title = 'Class',
                bbox_to_anchor = (0, 0), loc = 'best', 
                bbox_transform = plt.gcf().transFigure)

        handles = [Patch(facecolor = lut[name]) for name in lut]
        g.ax_col_dendrogram.legend(handles, lut, title = 'Tissues',
            bbox_to_anchor = (0, 1), loc = 'best', 
            bbox_transform = plt.gcf().transFigure)

        handlesY = [Patch(facecolor = lutY[name]) for name in lutY]
        g.ax_col_dendrogram.legend(handlesY, lutY, title = 'Samples',
            bbox_to_anchor = (0, 1), loc = 'best', 
            bbox_transform = plt.gcf().transFigure)

        g.savefig(str(i), dpi = 300)
        plt.close('all')
        c += 1

        # Free memory
        del(count)
        del(g)
        gc.collect()

if __name__ == '__main__':
    
    clustering_genes(
        cv_list = Config.args.list,
        g_corr = Config.g_corr,
        g_corr_rand = Config.args.corrRandG,
        g_corr_tissue = Config.args.corrTissueG,
        g_clustermaps = Config.g_clustermaps,        
        g_img_clstrRand = Config.args.IclstrRandG,
        g_img_clstrTissues = Config.args.IclstrTissueG)
    
    clustering_samples_genes(
        meta = Config.args.meta,
        cv_list = Config.args.list,
        counts = Config.counts,
        rand = Config.args.rand,
        by_tissue = Config.args.tissue,
        sg_clustermaps = Config.sg_clustermaps,
        sg_img_clstrRand = Config.args.IclstrRandSG,
        sg_img_clstrTissues = Config.args.IclstrTissueSG)
    
    clustering_samples(
        meta = Config.args.meta,
        s_corr = Config.s_corr,
		s_corr_rand = Config.args.corrRandS,
        s_corr_tissues = Config.args.corrTissueS,
        s_clustermaps = Config.s_clustermaps,
        s_img_clstrRand = Config.args.IclstrRandS,
        s_img_clstrTissues = Config.args.IclstrTissueS)