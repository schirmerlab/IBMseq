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




#load merged_object
merge_path="your_merge_path"
adata=sc.read_h5ad(merge_path)


# ## Create qc plots
:


#check sample annotation
adata.obs['sample_id'].cat.categories




#create dictionary with individual anndata objects per sample
sample_dic={}
for sample in adata.obs['sample_id'].cat.categories:
    sample_dic[sample]=adata[adata.obs['sample_id'].isin([sample])]

#create list with sample names
samples=[x for x in sample_dic]



#plot and save qc plot per sample
for sample in samples:
    sc.pl.scatter(sample_dic[sample], x='nCount_RNA' , y='nFeature_RNA', color='percent.mt',  save=f"{sample}_qc_pre_fil.pdf")


# ## Filter sampels 


#create dictionary with filtered samples 
qc_dic={}
for item in sample_dic:
    x=sample_dic[item]
    
    #choose cutoffs after checking qc plots
    x=x[x.obs['nCount_RNA']>1000]
    x=x[x.obs['nCount_RNA']<20000]
    x = x[x.obs['percent.mt'] < 2]
    x = x[x.obs['nFeature_RNA'] > 800]
    qc_dic[item]=x


# ## Doublet detection using scrublet 

import scrublet as scr



#predefine doublet detection function
#here we use a fixed threshold alt. automatic detection can be used
def pp_dt(adata):
    adata.var_names_make_unique()
    adata.raw=adata.copy()
    scrub= scr.Scrublet(adata.X , expected_doublet_rate=0.07)
    adata.obs['doublet_scores_scrublet'], adata.obs['predicted_doublets_scrublet'] = scrub.scrub_doublets(min_gene_variability_pctl=85, n_prin_comps=30)
    db_array=scrub.call_doublets(threshold=0.25)
    scrub.plot_histogram()
    adata.obs['predicted_doublets_scrublet']=db_array
    adata_unnorm=adata.raw.to_adata()
    adata_unnorm.obs=adata.obs
    return adata


#run doublet detection per sample
dt_dic={}
for smpl in qc_dic:
    dt_dic[smpl]= pp_dt(qc_dic[smpl])



#only include singlets
cleaned_dic={}
for it in dt_dic:
    cleaned_dic[it]= dt_dic[it][dt_dic[it].obs['predicted_doublets_scrublet']==False]


# ## Create qc plots after filtering  


#create qc plots after filtering if needed
samples=[x for x in cleaned_dic]
for sample in samples:
    sc.pl.scatter(cleaned_dic[sample], x='nCount_RNA' , y='nFeature_RNA', color='percent.mt',  save=f"{sample}_qc_post_fil.pdf")


# ### Create merged object 



#merged individual qc sampels together and filter genes
merged = sc.concat( cleaned_dic, join='outer')
sc.pp.filter_genes(merged, min_cells = 10)



#define path and save output#
merged_output_path='your_output_dictionary'
merged.write(merged_output_path)

