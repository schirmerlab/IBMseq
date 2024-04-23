#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import anndata


# ## Load Objects 



#load integrated spatial transcriptomics object
spatial_path='spatial_int_path'
adata=sc.read_h5ad(spatial_path)




#load deconvulution results for subsets
adata_immune=adata.copy()
adata_file = "immune_deconvolution_output"
adata_immune_met = sc.read_h5ad(f"{adata_file}/sp.h5ad")
adata_immune.obsm=adata_immune_met.obsm




adata_muscle=adata.copy()
adata_file = "muscle_deconvolution_output"
adata_muscle_met = sc.read_h5ad(f"{adata_file}/sp.h5ad")
adata_muscle.obsm=adata_muscle_met.obsm




adata_stroma=adata.copy()
adata_file = "stroma_deconvolution_output"
adata_stroma_met = sc.read_h5ad(f"{adata_file}/sp.h5ad")
adata_stroma.obsm=adata_stroma_met.obsm


# ## Compute correlation 


#create correlation function
def correlation_dic(adata):
    clusters=[x.split('sf_')[-1] for x in adata.obsm['q05_cell_abundance_w_sf'].columns]
    # get proportion mappings per celltype in a joint dic
    prop_dic=dict()
    for index , column in enumerate(adata.obsm['q05_cell_abundance_w_sf']):
        prop_dic[clusters[index]]=np.array(adata.obsm['q05_cell_abundance_w_sf'][column])
    return(prop_dic)


#run correlation
dic_muscle=correlation_dic(adata_muscle)
dic_immune=correlation_dic(adata_immune)
dic_stroma=correlation_dic(adata_stroma)




dic_immune={k : dic_immune[k] for k in clusters}



#merge Data Frames
dic_immune.update(dic_muscle)
dic_immune.update(dic_stroma)
prop_dic=dic_immune



#create dataframe with correlation values
#create count for x axis(column)
df_correlations=pd.DataFrame()
df_p_values=pd.DataFrame()
for x_index in prop_dic:
    value_x=prop_dic[x_index]
    y_values=[]
    y_pval=[]
    for y_index in prop_dic:
        value_y=prop_dic[y_index]
        y_values.append(stats.pearsonr(value_x, value_y)[0])
        y_pval.append(stats.pearsonr(value_x, value_y)[1])
    df_correlations[x_index]=y_values
    df_p_values[x_index]=y_pval
df_correlations['index']=prop_dic.keys()
df_correlations=df_correlations.set_index('index')


# ## Plot correlation  



#create plot with celltype corelation per spot fig 4d
fig, ax = plt.subplots()
plt.imshow(df_correlations, cmap='RdBu_r', interpolation='nearest', vmin=-1)
xlabel=prop_dic.keys()
ylabel= prop_dic.keys()
#plt.legend()
ax.set_xticklabels(xlabel)
plt.xticks(np.arange(len(prop_dic.keys())))
plt.yticks(np.arange(len(prop_dic.keys())))
ax.set_yticklabels(ylabel)
plt.colorbar()
plt.xticks(rotation=90)
#plt.savefig('/home/thomast/coding/cell2location_plots_2/correlation/per_spot_correlation_subtypes_scale.pdf', bbox_inches='tight', dpi=500)
plt.show()

