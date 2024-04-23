#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import bbknn
import anndata as ad
import os



#change figuresize
sc.set_figure_params(figsize=(8,6))


# ## Define functions 



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



#function for bbknn integration
def int_bbknn(adata, n_pcs=25, batch_ident='sample_id', ridge_reg=True ):
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


# ## Loading and normalization 



#load qc object
qc_path="your_qc_object_path"
merged=sc.read_h5ad(qc_path)


#log transfrom and counts storing
merged.layers['counts']=merged.X.copy()
sc.pp.normalize_total(merged, target_sum=1e4)
sc.pp.log1p(merged)
merged.raw=merged



merged=calc_hvg(merged,)




#regression of techical covariats if needed and scaling
sc.pp.regress_out(merged,['nCount_RNA', 'percent.mt'])
sc.pp.scale(merged, max_value=10)


# ## Batch correction using bbknn



merged=int_bbknn(merged)



#check integratin results
sc.pl.umap(merged, color=['sample_id','condition','batch'])




sc.tl.leiden(merged, resolution=0,6)




sc.pl.umap(merged, color='leiden')



#scan clusters markers and display top hits
sc.tl.rank_genes_groups(merged, groupby='leiden', key_added= "rank_genes_groups")
pd.DataFrame(merged.uns['rank_genes_groups']['names']).head(20)



#rename clusters based on markergenes here exp. overview clustering
cluster_names=['Type 2 MN','FAPs','Type 1 MN','Reactive MN','MФ','T cells','Satellite cells','Endothelial cells','Pericytes','DC','Damaged MN','Adipocytes','naFb','Muscle spindle']
merged.rename_categories('leiden', cluster_names)


# ## Create dotplots



#reorder clusters based on biological groups if needed 
order=['Type 1 MN','Type 2 MN','Reactive MN','Damaged MN','Satellite cells','Muscle spindle','FAPs','naFb','Endothelial cells','Pericytes','Adipocytes','T cells','DC','MФ']
merged.obs["leiden_ordered"] = (
    merged.obs["leiden"]
    .cat.reorder_categories(order
    )
)



#create dotplot
genes=['MRC1','CADM1','THEMIS','ADIPOQ','RGS5','VWF','TENM2','PDGFRA','PIEZO2','PAX7','GADD45A','COL19A1','ATP2A1','ATP2A2']
sc.pl.dotplot(merged, genes, 'leiden_ordered', save=True)

