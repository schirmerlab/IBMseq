setwd("/your/working/directory")

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(ggh4x)

# R version 4.2.2

######## Differential abundance analysis on different subsets using miloR. For more details on how to use this package please see "https://github.com/MarioniLab/miloR" ########

#load your files

subset_lognorm_var_genes <- readRDS("your/path/your_subset_lognorm_var_genes.rds")
subset_lognorm_all_genes <- readRDS("your/path/your_subset_lognorm_all_genes.rds")

subset_lognorm_all_genes <- logNormCounts(subset_lognorm_all_genes, assay.type = "X")

reducedDim(subset_lognorm_all_genes, "PCA") <- reducedDim(subset_lognorm_var_genes, "PCA")
reducedDim(subset_lognorm_all_genes, "UMAP") <- reducedDim(subset_lognorm_var_genes, "UMAP")

#create miloR object

milo_subset <- Milo(subset_lognorm_all_genes)


set.seed(8)

milo_subset <- buildGraph(milo_subset, k = 20, d = 20)
milo_subset <- makeNhoods(milo_subset, prop = 0.05, k = 20, d=20, refined = TRUE)

plotNhoodSizeHist(milo_subset)
milo_subset <- countCells(milo_subset, meta.data = data.frame(colData(milo_subset)), samples="sample_id")
head(nhoodCounts(milo_subset))

#set up a model and define contrasts

subset_design <- data.frame(colData(milo_subset))[,c("sample_id", "condition")]
subset_design <- distinct(subset_design)
rownames(subset_design) <- subset_design$sample_id
subset_design <- subset_design[colnames(nhoodCounts(milo_subset)), , drop=FALSE]
table(subset_design$condition)

model <- model.matrix(~ 0 + condition, data=subset_design)
contrast.IBMvCTRL <- c("conditionIBM - conditionCTRL")
contrast.IMNMvCTRL <- c("conditionIMNM - conditionCTRL")
contrast.IBMvIMNM <- c("conditionIBM - conditionIMNM") 

#### IBM v CTRL ####
mod.constrast.IBMvCTRL <- makeContrasts(contrasts=contrast.IBMvCTRL, levels=model)

mod.constrast.IBMvCTRL
set.seed(8)

da_results_subset_IBMvCTRL <- testNhoods(milo_subset, design = ~ 0 + condition, design.df = subset_design, model.contrasts = contrast.IBMvCTRL, fdr.weighting="graph-overlap")

table(da_results_subset_IBMvCTRL$SpatialFDR < 0.1)


milo_subset <- buildNhoodGraph(milo_subset)
plotUMAP(milo_subset, color_by="leiden_ordered", point_size = 0.1) + plotNhoodGraphDA(milo_subset, da_results_subset_IBMvCTRL, alpha=0.05) +
  plot_layout(guides="keep")

ggplot(da_results_subset_IBMvCTRL, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)

da_results_subset_IBMvCTRL <- annotateNhoods(milo_subset, da_results_subset_IBMvCTRL, coldata_col = "leiden_ordered")
head(da_results_subset_IBMvCTRL)
ggplot(da_results_subset_IBMvCTRL, aes(leiden_ordered_fraction)) + geom_histogram(bins=50)
da_results_subset_IBMvCTRL$leiden_ordered <- ifelse(da_results_subset_IBMvCTRL$leiden_ordered_fraction < 0.7, "Mixed", da_results_subset_IBMvCTRL$leiden_ordered)

plotDAbeeswarm(da_results_subset_IBMvCTRL, group.by = "leiden_ordered", )

#### IBM v IMNM #### 

mod.contrast.IBMvIMNM <- makeContrasts(contrasts=contrast.IBMvIMNM, levels=model)


mod.contrast.IBMvIMNM
set.seed(8)

da_results_subset_IBMvIMNM <- testNhoods(milo_subset, design = ~ 0 + condition, design.df = subset_design, model.contrasts = contrast.IBMvIMNM, fdr.weighting="graph-overlap")


table(da_results_subset_IBMvIMNM$SpatialFDR < 0.1)


milo_subset <- buildNhoodGraph(milo_subset)
plotUMAP(milo_subset, color_by="leiden_ordered", point_size = 0.1) + plotNhoodGraphDA(milo_subset, da_results_subset_IBMvIMNM, alpha=0.05) +
  plot_layout(guides="keep")
ggplot(da_results_subset_IBMvIMNM, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)

da_results_subset_IBMvIMNM <- annotateNhoods(milo_subset, da_results_subset_IBMvIMNM, coldata_col = "leiden_ordered")
head(da_results_subset_IBMvIMNM)
ggplot(da_results_subset_IBMvIMNM, aes(leiden_ordered_fraction)) + geom_histogram(bins=50)
da_results_subset_IBMvIMNM$leiden_ordered <- ifelse(da_results_subset_IBMvIMNM$leiden_ordered_fraction < 0.7, "Mixed", da_results_subset_IBMvIMNM$leiden_ordered)

plotDAbeeswarm(da_results_subset_IBMvIMNM, group.by = "leiden_ordered")

#### IMNM v CTRL #### 

mod.contrast.IMNMvCTRL <- makeContrasts(contrasts=contrast.IMNMvCTRL, levels=model)

mod.contrast.IMNMvCTRL
set.seed(8)

da_results_subset_IMNMvCTRL <- testNhoods(milo_subset, design = ~ 0 + condition, design.df = subset_design, model.contrasts = contrast.IMNMvCTRL, fdr.weighting="graph-overlap")

table(da_results_subset_IMNMvCTRL$SpatialFDR < 0.1)

milo_subset <- buildNhoodGraph(milo_subset)
plotUMAP(milo_subset, color_by="leiden_ordered", point_size = 0.1) + plotNhoodGraphDA(milo_subset, da_results_subset_IMNMvCTRL, alpha=0.05) +
  plot_layout(guides="keep")

ggplot(da_results_subset_IMNMvCTRL, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)



da_results_subset_IMNMvCTRL <- annotateNhoods(milo_subset, da_results_subset_IMNMvCTRL, coldata_col = "leiden_ordered")
head(da_results_subset_IMNMvCTRL)
ggplot(da_results_subset_IMNMvCTRL, aes(leiden_ordered_fraction)) + geom_histogram(bins=50)
da_results_subset_IMNMvCTRL$leiden_ordered <- ifelse(da_results_subset_IMNMvCTRL$leiden_ordered_fraction < 0.7, "Mixed", da_results_subset_IMNMvCTRL$leiden_ordered)

plotDAbeeswarm(da_results_subset_IMNMvCTRL, group.by = "leiden_ordered")


#### multiple comparison #### 

#create empty data frame

da_results_subset_IBMvCTRL_adj <- data.frame(matrix(ncol = 9, nrow = nrow(da_results_subset_IBMvCTRL)))
colnames(da_results_subset_IBMvCTRL_adj) <- colnames(da_results_subset_IBMvCTRL)
da_results_subset_IBMvCTRL_adj[, "SpatialFDR_old"] = NA

da_results_subset_IBMvIMNM_adj <- data.frame(matrix(ncol = 9, nrow = nrow(da_results_subset_IBMvIMNM)))
colnames(da_results_subset_IBMvIMNM_adj) <- colnames(da_results_subset_IBMvIMNM)
da_results_subset_IBMvIMNM_adj[, "SpatialFDR_old"] = NA

da_results_subset_IMNMvCTRL_adj <- data.frame(matrix(ncol = 9, nrow = nrow(da_results_subset_IMNMvCTRL)))
colnames(da_results_subset_IMNMvCTRL_adj) <- colnames(da_results_subset_IMNMvCTRL)
da_results_subset_IMNMvCTRL_adj[, "SpatialFDR_old"] = NA

#adjust Spatial FDR values as there are three conditions
#adjusted Spatial FDR values will be in a column called "SpatialFDR" and non-adjusted values will be in the "SpatialFDR_old" column. 

for (i in 1:nrow(da_results_subset_IBMvCTRL)) {
  
  nam <- paste("df_IBMvCTRL", i, sep = "")
  nam2 <- paste("df_IBMvIMNM", i, sep = "")
  nam3 <- paste("df_IMNMvCTRL", i, sep = "")
  nam4 <- paste("da_results_subset_mc", i, sep = "")
  assign(nam, da_results_subset_IBMvCTRL[da_results_subset_IBMvCTRL$Nhood==i,])
  assign(nam2, da_results_subset_IBMvIMNM[da_results_subset_IBMvIMNM$Nhood==i,])
  assign(nam3, da_results_subset_IMNMvCTRL[da_results_subset_IMNMvCTRL$Nhood==i,])
  
  df <- rbind(get(nam), get(nam2), get(nam3))
  df$SpatialFDR_old <- df$SpatialFDR
  df$SpatialFDR <- p.adjust(df$SpatialFDR, method = "BH")
  da_results_subset_IBMvCTRL_adj[i,] <- df[1,]
  da_results_subset_IBMvIMNM_adj[i,] <- df[2,]
  da_results_subset_IMNMvCTRL_adj[i,] <- df[3,]
  rm(list = ls(pattern = "^df"))
  rm(list = ls(pattern = "^nam"))
  
}

#If you want, you can check how many neighbourhoods remained significant

table(da_results_subset_IBMvCTRL_adj$SpatialFDR < 0.1)
table(da_results_subset_IBMvIMNM_adj$SpatialFDR < 0.1)
table(da_results_subset_IMNMvCTRL_adj$SpatialFDR < 0.1)

#arrange cell types for plot generation

da_results_subset_IBMvCTRL_adj$leiden_ordered <- factor(da_results_subset_IBMvCTRL_adj$leiden_ordered, levels = c("Mixed", "subset_cluster1", "etc."))
levels(da_results_subset_IBMvCTRL_adj$leiden_ordered)
table(factor(da_results_subset_IBMvCTRL_adj$leiden_ordered))

da_results_subset_IBMvIMNM_adj$leiden_ordered <- factor(da_results_subset_IBMvIMNM_adj$leiden_ordered, levels = c("Mixed", "subset_cluster1", "etc."))
levels(da_results_subset_IBMvIMNM_adj$leiden_ordered)
table(factor(da_results_subset_IBMvIMNM_adj$leiden_ordered))

da_results_subset_IMNMvCTRL_adj$leiden_ordered <- factor(da_results_subset_IMNMvCTRL_adj$leiden_ordered, levels = c("Mixed", "subset_cluster1", "etc."))
levels(da_results_subset_IMNMvCTRL_adj$leiden_ordered)
table(factor(da_results_subset_IMNMvCTRL_adj$leiden_ordered))

#plots

plotDAbeeswarm(da_results_subset_IBMvCTRL_adj, group.by = "leiden_ordered") + force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))
plotDAbeeswarm(da_results_subset_IBMvIMNM_adj, group.by = "leiden_ordered") + force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))
plotDAbeeswarm(da_results_subset_IMNMvCTRL_adj, group.by = "leiden_ordered") + force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))

plotNhoodGraphDA(milo_subset, da_results_subset_IBMvCTRL_adj, alpha=0.1) + force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))
plotNhoodGraphDA(milo_subset, da_results_subset_IBMvIMNM_adj, alpha=0.1) + force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))
plotNhoodGraphDA(milo_subset, da_results_subset_IMNMvCTRL_adj, alpha=0.1) + force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))
  
plotDAbeeswarm(da_results_subset_IBMvCTRL_adj, group.by = "leiden_ordered") + ggtitle("IBM - CTRL") +
  plotDAbeeswarm(da_results_subset_IBMvIMNM_adj, group.by = "leiden_ordered") + ggtitle("IBM - IMNM") +
  plotDAbeeswarm(da_results_subset_IMNMvCTRL_adj, group.by = "leiden_ordered") + ggtitle("IMNM - CTRL") +
  plotNhoodGraphDA(milo_subset, da_results_subset_IBMvCTRL_adj, alpha=0.1) +
  plotNhoodGraphDA(milo_subset, da_results_subset_IBMvIMNM_adj, alpha=0.1) +
  plotNhoodGraphDA(milo_subset, da_results_subset_IMNMvCTRL_adj, alpha=0.1) +
  plotUMAP(milo_subset, color_by="leiden_ordered", point_size = 0.1) +
 plot_layout(guides="keep", ncol = 3)

#### save results #### 

saveRDS(da_results_subset_IBMvCTRL_adj, "da_results_subset_IBMvCTRL_adj.rds")
saveRDS(da_results_subset_IBMvIMNM_adj, "da_results_subset_IBMvIMNM_adj.rds")
saveRDS(da_results_subset_IMNMvCTRL_adj, "da_results_subset_IMNMvCTRL_adj.rds")
saveRDS(milo_subset, "milo_subset.rds")

