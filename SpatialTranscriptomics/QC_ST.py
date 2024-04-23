#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import decoupler as dc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import bbknn
import anndata as ad
import scvi
import os
import squidpy as sq



#change figuresize
sc.set_figure_params(figsize=(8,6))




#load individual h5ad
parent_path='your_path_to_h5ads'
sample_names=os.listdir(parent_path)


# ## Sample reading 



#create local objects if needed
for sample in sample_names:
    locals()[sample[:]]=sc.read_h5ad(f"{parent_path}/{sample}")


# ## Add meta information



#create dictionary with each sample
sample_dic=dict()
for sample in sample_id:
    sample_dic[sample]=locals()[sample]
#sample_id_s=sorted(sample_id)




#create qc plots per sample
for sample in sample_dic:
    sample_dic[sample].var["mt"] = sample_dic[sample].var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(sample_dic[sample], qc_vars=["mt"], inplace=True)
    sc.pl.scatter(sample_dic[sample], x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', save=f"{sample}_qc_pre_fil.pdf")



#Check additional qc metrics if needed
qc_dic={}
for item in sample_id_s[:]:
    x=sample_dic[item]
    
    x.var["mt"] = x.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(x, qc_vars=["mt"], inplace=True)
    #print('qc plots for ', i)
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.distplot(x.obs["total_counts"], kde=False, ax=axs[0])
    sns.distplot(x.obs["total_counts"][x.obs["total_counts"] < 4000], kde=False, bins=40, ax=axs[1])
    sns.distplot(x.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    sns.distplot(x.obs["n_genes_by_counts"][x.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
    qc_dic[item]=x


# # Filter samples


##create dictionary with filtered samples 
qc_dic={}
for item in sample_dic:
    print(item)
    x=sample_dic[item]
    strt=x.n_obs
    print(f"#cells before filter: {x.n_obs}")
    sc.pp.filter_cells(x, min_counts=600)
    sc.pp.filter_cells(x, min_genes=200)
    print(f"#cells after MT filter: {x.n_obs}")
    sc.pp.filter_genes(x, min_cells=10)
    end=x.n_obs
    print('pct change ', end/strt)
    qc_dic[item]=x


# ## Create qc plots after filtering  



#additional qc pltos after filtering if needed
for sample in qc_dic:
    sc.pl.scatter(qc_dic[sample], x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', save=f"{sample}_qc_post_fil.pdf")


# # Save results 



#make indexes unique and save individual samples
p_path='your_output_path'
for sample in sample_id:
    x=sample_dic[sample]
    x.var_names_make_unique()
    x.obs_names_make_unique()
    x.write(p_path + sample + '_qc.h5ad')


# # Merge object 



#create merged object and save object
merged=sc.concat(qc_dic, join='outer', label='sample_id', uns_merge='unique',
    index_unique="-")
merged.write(p_path + 'merged' + '_qc.h5ad')

