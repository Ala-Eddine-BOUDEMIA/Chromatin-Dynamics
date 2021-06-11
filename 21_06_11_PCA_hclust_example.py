#!/bin/env python3

import numpy as np
import scipy
import pandas as pd
import matplotlib
import seaborn as sns
import sklearn

#####################################################################################################################
#PCA plot

from sklearn.decomposition import PCA

#SAMPLES TO CMP PCA AND HCLUST
samples_for_pca = samples_lr_z_10kb[qc_list]
#samples are organised in a table where every row is a genomic bin and every column is a sample

sample_details_for_pca = sample_IP_H3vars_filt.loc[qc_list]
#sample details are organised in a table where every row is a sample and every column is a sample attribute (e.g. input/IP, read number, antibody used, sample condition, etc.)

pca = PCA().fit(samples_for_pca.dropna().T)
Y = pca.fit_transform(samples_for_pca.dropna().T)
ratio = pca.explained_variance_ratio_ * 100

#plt.scatter(Y[:, 0], Y[:, 1])

d = pd.DataFrame(index = samples_for_pca.T.index) #set column names (samples) of the sample df as index of d
d["PC1"] = Y[:, 0]
d["PC2"] = Y[:, 1]
d = d.join(sample_details_for_pca)
#add sample details (every row is a sample, every col is an attribute)


for hue in ['IP', 'antibody', 'study_accession']:
    sns.scatterplot(x = "PC1", y = "PC2", hue = hue, data = d) #color points on the plot by the respective attribute
    plt.xlabel('PC1 %s' %(round(ratio[0],2)))
    plt.ylabel('PC2 %s' %(round(ratio[1],2)))
    plt.savefig('%s/processed/mouse/DNA/enrichment_per_chrs/PCA_of_steady_state_H3vars_hue_%s.png' %(basefolder, hue), dpi = 300)
    plt.close('all')
    gc.collect()
#[plt.text(x = x, y = y, s = s) for x, y, s in d[["PC1", "PC2", "sample_title"]].to_numpy()]

#Hierachical clustering w/ Pearson/Spearman R
#for QC list of samples
f = sns.clustermap(samples_lr_z_10kb[qc_list].corr(), cmap = "YlGnBu", vmin = 0, vmax = 1,
                   xticklabels = sample_IP_H3vars_filt.loc[qc_list]['sample_title'],
                   yticklabels = sample_IP_H3vars_filt.loc[qc_list]['sample_title'],
                   #row_colors = [red if p == 'FLAG' else yellow if p == 'HA' else blue if p == 'SNAP' else dark_green if p == 'GFP' else
                    #             dark_purple if p == 'Millipore 09-838' else greyish if p == 'Abcam ab92628' else blueish if p == 'Abcam ab1791' else pale_yellow if p == 'Millipore ABE154' else light_green #if p == 'Active Motif 91191'
                    #             for p in sample_IP_H3vars_filt.loc[qc_list]['antibody']],
                   row_colors = [red if p == 'PRJDB2002' else yellow if p == 'PRJNA117667' else blue if p == 'PRJNA254837' else dark_green if p == 'PRJNA314327' else
                                dark_purple if p == 'PRJNA254708' else greyish if p == 'PRJNA627246' else blueish if p == 'PRJNA580180' else pale_yellow if p == 'PRJDB3075' else light_green #if p == 'PRJNA471722'
                                for p in sample_IP_H3vars_filt.loc[qc_list]['study_accession']],
                   col_colors = [green if s == 'H3.3' else light_green if s == 'H3.1S31p' else purple if s == 'H3.1' else light_purple if (s == 'H3.2') | (s == 'H3.1S31') | (s == 'H3.1/2') else blue #if s == 'H3gen'
                                 for s in sample_IP_H3vars_filt.loc[qc_list]['IP']],
                   figsize = [25, 25], fmt = '.2f', annot = True, method = 'complete')

plt.savefig('%s/processed/mouse/DNA/enrichment_per_chrs/hclust_of_all_H3vars_qc_ChIP_rows_prj_cols_H3var_annot.png' %(basefolder), dpi = 300)
plt.close('all')
gc.collect()
