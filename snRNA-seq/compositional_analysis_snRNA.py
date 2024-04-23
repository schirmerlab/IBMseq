#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
pd.set_option('display.max_columns', 50)




import pickle as pkl
import sccoda
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import anndata as ad



#change figuresize
sc.set_figure_params(figsize=(8,6))




#load merged integrated object
integrated_path='your_integrated_path'
merged=sc.read_h5ad(integrated_path)



#plot umap overview if needed
sc.pl.umap(merged, color='leiden_ordered')


#create metadata data Frame for sccoda loading, choose categories of interest
data_per_sample=merged.obs.drop_duplicates(subset='sample_id') 
meta_data=data_per_sample[['sample_id','sex','age','condition']]
cov_df=meta_data
cov_df['index']=cov_df['sample_id']
cov_df=cov_df.set_index('index')
cov_df['sample_n']=cov_df['sample_id']
del(cov_df['sample_id'])
newlist=[str(x) for x in cov_df['condition']]
cov_df['condition']=newlist



#create sccoda object
data_scanpy_1=dat.from_scanpy(merged ,cell_type_identifier='leiden_ordered', sample_identifier='sample_id', covariate_df=cov_df)


# ## Compositional analysis 



#create default plots depending on question
# Stacked barplot for each sample
viz.stacked_barplot(data_scanpy_1, feature_name="sample_n")

#plt.savefig(your_def_path.pdf', dpi=500, format='pdf')
plt.legend( loc='lower right', fontsize='xx-small')
plt.show()

# Stacked barplot for the levels of "Condition", plot used in fig. 1
viz.stacked_barplot(data_scanpy_1, feature_name="condition")
plt.legend( loc='lower right', fontsize='xx-small')

#plt.savefig('your_def_path.pdf', dpi=500, format='pdf')
plt.show()



#create additional boxplots per cluster if needed
#choose between y_scale "log" or "relative"
viz.boxplots(
    data_scanpy_1,
    feature_name="condition",
    plot_facets=False,
    y_scale="log",
    add_dots=False,
    cmap=['lightblue',  'orange','red',],
    level_order=['CTRL','IMNM','IBM'],
)
plt.savefig('your_compositional_plot_path', dpi=500, format='pdf')
plt.show()

