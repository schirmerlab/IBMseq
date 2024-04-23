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




#create abundance function
def abundance_pipeline(adata, ):
    dic_cluster_barcodes=dict()
    for cluster in adata.obs.leiden.cat.categories:
        barcodes=[x for x in adata.obs.leiden[adata.obs.leiden.values == str(cluster)].index]
        dic_cluster_barcodes[cluster]=barcodes

    dic_df=dict()
    for key in dic_cluster_barcodes:
        df=adata.obsm['q05_cell_abundance_w_sf'][adata.obsm['q05_cell_abundance_w_sf'].index.isin(dic_cluster_barcodes[key])]
        dic_df[key]=df

    dic_mean_per_cluster=dict()
    for df in dic_df:
        temp_dic=dict()
        columns=list(dic_df[df].columns)
        clusters=[x.split('sf_')[-1] for x in columns]
        for index, cluster in enumerate(clusters):
            temp_dic[cluster]=dic_df[df][columns[index]].mean()
        dic_mean_per_cluster[df]=temp_dic

    df_proportions=pd.DataFrame()

    for key in dic_mean_per_cluster:
        #count=0
        df_proportions=df_proportions.append(dic_mean_per_cluster[key] ,ignore_index=True)
    return(df_proportions)




#calculate abundance per subset and combine Data Frames
df_proportions_immune=abundance_pipeline(adata_immune)
df_proportions_muscle=abundance_pipeline(adata_muscle)
df_proportions_stroma=abundance_pipeline(adata_stroma)

df_combined=df_proportions_immune.join(df_proportions_muscle)
df_proportions=df_combined.join(df_proportions_stroma)


# ## Plotting values


#craete raw plot
fig, ax = plt.subplots()
plt.imshow(df_proportions, cmap='viridis', interpolation='nearest')
xlabel=(list(df_proportions.columns))
ylabel= list(df_proportions.index)
#plt.legend()
ax.set_xticklabels(xlabel)
plt.xticks(np.arange(27))
plt.yticks(np.arange(5))
ax.set_yticklabels(ylabel)
plt.xticks(rotation=90)
plt.show()


# ## Standardized values



from sklearn.preprocessing import StandardScaler


#scale values using sklearn StandardScaler
std_scaler = StandardScaler() 
df_scaled = std_scaler.fit_transform(df_proportions.to_numpy())
df_scaled = pd.DataFrame(df_scaled, columns=df_proportions.columns)


# # Spot proportions scaled 



#fig 4d heatmap with scaled abundances
fig, ax = plt.subplots()
plt.imshow(df_scaled, cmap='RdBu_r', interpolation='nearest', vmin=-2, vmax=2)
xlabel=(list(df_scaled.columns))
ylabel= list(df_scaled.index)
ax.set_xticklabels(xlabel)
plt.xticks(np.arange(27))
plt.yticks(np.arange(5))
plt.colorbar()
ax.set_yticklabels(ylabel)
plt.xticks(rotation=90)
#plt.savefig('heatmap_abundances_path.pdf', bbox_inches='tight', dpi=500)
plt.show()

