#!/usr/bin/env python
# coding: utf-8



import scanpy as sc
import anndata
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import bbknn
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42


# ## Load objects



#load merged qc object
merged_path='your_merged_spatial_path'
adata_vis_unnorm=sc.read_h5ad(merged_path)




#load cell2location results
adata_file = 'your_deconvolution_results_path'
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load('your_deconvolution_output', adata_vis)




#formatting for compatability
adata_vis.uns['spatial']=adata_vis_unnorm.uns['spatial']
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
adata_vis.obs['sample']=[ s.upper() for s in adata_vis.obs['sample']]


# ## Create custom functions 



#function for normalization and log transformation
def norm_log(adata, layers='counts'):
    adata.layers[layers]=adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw=adata
    return(adata)



#function for calculating highly variable genes taking batches into acount
def calc_hvg(adata, n_hvg= 2500, inplace=True,  batch_key='sample_id'):
    sc.pp.highly_variable_genes(adata, batch_key=batch_key)
    adata.var=adata.var[['highly_variable','highly_variable_nbatches']]
    batch_msk=np.array(adata.var.highly_variable_nbatches >1)
    hvg=adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(n_hvg).index
    adata.var['highly_variable']=[g in hvg for g in adata.var.index]
    adata.var=adata.var[['highly_variable','highly_variable_nbatches']]
    if inplace == True:
        adata=adata[:,hvg]
    return(adata)



#function for scaling and regression if needed
def reg_scale(adata, reg=True, reg_parm=['total_counts', 'pct_counts_mt'], scale_max=10):
    if reg is True:
        sc.pp.regress_out(adata, reg_parm)
    sc.pp.scale(adata, max_value=scale_max)
    return(adata)



#function for bbknn integration
def int_bbknn(adata, n_pcs=15, batch_ident='sample_id', ridge_reg=True ):
    sc.pp.pca(adata, )
    sc.pl.pca_variance_ratio(adata)
    plt.show()
    bbknn.bbknn(adata, batch_key=batch_ident,n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=batch_ident, )
    plt.show()
    sc.tl.leiden(adata, resolution=0.6)
    sc.pl.umap(adata, color = 'leiden')
    if ridge_reg is True:
        print('running regression')
        bbknn.ridge_regression(adata, batch_key=batch_ident, confounder_key=['leiden'])
        sc.pp.pca(adata)
        bbknn.bbknn(adata, batch_key=batch_ident, n_pcs=n_pcs)
        sc.tl.umap(adata)
    sc.pl.umap(adata, color=batch_ident)
    plt.show()
    return(adata)


# ## Cell2location per slide plotting



from cell2location.utils import select_slide
#load each slide
slide_c56 = select_slide(adata_vis, 'C56')
slide_c50 = select_slide(adata_vis, 'C50')
slide_c47 = select_slide(adata_vis, 'C47')
slide_ibm29 = select_slide(adata_vis, 'IBM29')
slide_ibm31 = select_slide(adata_vis, 'IBM31')
slide_ibm35 = select_slide(adata_vis, 'IBM35')
slide_IMNM1=select_slide(adata_vis, 'IMNM1')
slide_IMNM4=select_slide(adata_vis, 'IMNM4')


#viszalize deconvolution results on sldie fig 4e
#here just one example
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(slide_c47, cmap='magma',
                  color=['Type 2 MN'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  vmin=0, vmax='p99.2'
                 )


# ## Joint analysis


#run normalization, hvg calculation 
adata_vis_unnorm=adata_vis.copy()
adata_vis=norm_log(adata_vis)
adata_vis=calc_hvg(adata_vis,)
adata=reg_scale(adata_vis, reg=False)


#run bbk integration, choose pcs based on complexity of the dataset and elbow plot
adata_vis=int_bbknn(adata_vis, n_pcs=10)




#rerun clustering and check integration results
sc.tl.leiden(adata_vis , resolution=0.5)
sc.pl.umap(adata_vis, color=['condition', 'sample_id',])


#renaming clusters
cluster_names=['Niche 0','Niche 1','Niche 2', 'Niche 3', 'Niche 4']
data_vis.rename_categories('leiden', cluster_names)


# ## Compositional Table for composition
:


#calculate compositional table with samples as columns and clusters as rows
df_composition=pd.DataFrame()
for sample in adata_vis.obs.sample_id.cat.categories:
    subset=adata_vis[adata_vis.obs.sample_id.isin([sample])]
    value_index=list(subset.obs.leiden.value_counts().index)
    print(value_index)
    value=list(subset.obs.leiden.value_counts())
    print(value)
    sorted_values=[value for _, value in sorted(zip(value_index, value))]
    df_composition[sample]=sorted_values



#display table and save as csv
df_composition
#df_composition.to_csv('composition_visium_clusters.csv')


# ## Dotplot per cluster 



#create dotplot of genes of interes fig 5d
genes=['GADD45A','RNF7','NORAD','ACHE','CD3E','CXCL9']
sc.pl.stacked_violin(adata_vis, genes, groupby='leiden', swap_axes=True,save=True )


# ## Plotting ACHE expression on umap 



sc.pl.umap(adata_vis, color='ACHE', )


# ## Spatial neighbors heatmap


import squidpy as sq
#calculate spatial neighbors and enrichment
sq.gr.spatial_neighbors(adata_vis)
sq.gr.nhood_enrichment(adata_vis, cluster_key="leiden")
n_enrichment=adata_vis.uns['leiden_nhood_enrichment']['zscore']
pd.DataFrame(n_enrichment)
df=pd.DataFrame(n_enrichment)


#create spatial neighbors heatmap fig 4c
fig, ax = plt.subplots()
plt.imshow(n_enrichment, cmap='RdBu_r', interpolation='nearest', vmin=-12.5)
xlabel=adata_vis.obs.leiden.cat.categories
ylabel=adata_vis.obs.leiden.cat.categories
ax.set_xticklabels(xlabel)
ax.set_yticklabels(ylabel)
plt.xticks(np.arange(len(adata_vis.obs.leiden.cat.categories)))
plt.yticks(np.arange(len(adata_vis.obs.leiden.cat.categories)))
plt.colorbar()
plt.legend()
plt.xticks(rotation=90)
#plt.savefig('your_path_heatmap')
plt.show()



adta_vis.write('spatial_int_path')




