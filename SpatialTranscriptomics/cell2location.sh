#!/bin/env python
#SBATCH --partition=single
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=60gb
#SBATCH --gres=gpu:1



import scanpy as sc 
import anndata
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'
import cell2location
import scvi



#define paths 
visium_ref_path='path_int_ST_object'
single_cell_int_csv_path='path_csv_metadata_clustering_snRNA'
single_cell_unnorm_path='path_unnorm_snRNA_object'

#define output paths
basis='output_path'
regression_path=f"{basis}/regression_output"
deconvolution_path=f"{basis}/deconvolution_output"
qc_plot_path=f"{basis}/qc_plots"


#define parameters
cluster_key='leiden_ordered'
batch_key='sample_id'

#define hyperparameters
dt_alpha=200
n_cells_per_location=5


#load objects
adata_vis_unnorm=sc.read_h5ad(visium_ref_path)
adata_ref_unnorm=sc.read_h5ad(single_cell_unnorm_path)
df_cluster_info=pd.read_csv(single_cell_int_csv_path)


#part I run regression on cell abundance
#process visium
adata_vis_unnorm.obs['sample']=adata_vis_unnorm.obs['sample_id']
adata_vis_unnorm.obs_names_make_unique()
adata_vis_unnorm.var_names_make_unique()
adata_vis_unnorm.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis_unnorm.var.index]
adata_vis_unnorm.obsm['MT'] = adata_vis_unnorm[:, adata_vis_unnorm.var['MT_gene'].values].X.toarray()
adata_vis_unnorm = adata_vis_unnorm[:, ~adata_vis_unnorm.var['MT_gene'].values]


#process single cell
adata_ref_unnorm.obs=df_cluster_info
adata_ref=adata_ref_unnorm
adata_ref.var_names_make_unique()
adata_ref.var_names_make_unique()

#feature selection
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

#setup regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key=batch_key,
                        # cell type, covariate used for constructing signatures
                        labels_key=cluster_key,
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
			#categorical_covariate_keys=['Method']
		)

  
#create the regression model
mod = cell2location.models.RegressionModel(adata_ref)
mod.view_anndata_setup()

#run training
mod.train(max_epochs=250, use_gpu=True)

#export cell abundance
data_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)


# Save model
mod.save(regression_path, overwrite=True)
adata_file = f"{regression_path}/sc.h5ad"
#adata_ref.write(adata_file)

#try plotting qc
mod.plot_QC()
plt.savefig(f"{qc_plot_path}/regression")


#part II run cell2location deconvolution

#prepare data for cell2location
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']
adata_vis=adata_vis_unnorm



#find intersections in genes selected
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

#create decon model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=n_cells_per_location, 
    detection_alpha=dt_alpha
)
mod.view_anndata_setup()


#train deconv model
mod.train(max_epochs=30000,
          batch_size=None,
          train_size=1,
          use_gpu=True)
#export training plot
mod.plot_history(1000)
plt.savefig(f"{qc_plot_path}/decon_training")



#exporting results
# In this section, we export the estimated cell abundance (summary of the posterior distribution)
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(deconvolution_path, overwrite=True)
adata_file = f"{deconvolution_path}/sp.h5ad"
adata_vis.write(adata_file)

#export decon qcs
mod.plot_QC()
plt.savefig(f"{qc_plot_path}/decon_qc")
