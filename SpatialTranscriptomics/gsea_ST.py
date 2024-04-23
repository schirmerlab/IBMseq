#!/usr/bin/env python
# coding: utf-8


import gseapy as gp
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pld




#define inputs
single_cell_path= 'your_single_cell_object_path'
cluster_interest='your_cluster_of_interes' # exp. Type 2 MN
output_path='your_gsea_output_path'



adata=sc.read(single_cell_path)
adata.uns['log1p']['base']=None
adata=adata.raw.to_adata()




#define filter function for geneset
def filter_dictionary(dictionary, allowed_values):
    dictionary2=dict()
    for key, value in dictionary.items():
        dictionary2[key] = [item for item in value if item in allowed_values]
    return dictionary2




#create custom geneset with background considered
go_mdf=gp.get_library('GO_Biological_Process_2021',organism='Human')
background=[x for x in adata.var.index]

mod_gs=filter_dictionary(go_mdf, background




#subset adata object
bdata=adata.copy()
bdata.obs['clinterest']=bdata.obs['leiden_ordered']==cluster_interest
bdata.obs['clinterest'] = pd.Categorical(bdata.obs['clinterest'], categories=[True, False],ordered=True)
print('starting')
res = gp.gsea(data=bdata.to_df().T, # row -> genes, column-> samples
        gene_sets=mod_gs,
        cls=bdata.obs.clinterest,
        permutation_num=1000,
        permutation_type='phenotype',
        outdir=None, 
        method='s2n', # signal_to_noise
        threads= 16)
res.res2d.to_csv(f"{output_path}/{cluster_interest}_gsea.csv")
ptint('done')

