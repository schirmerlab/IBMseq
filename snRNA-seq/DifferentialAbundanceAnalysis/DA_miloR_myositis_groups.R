setwd("/your/working/directory")

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(tidyverse)
library(patchwork)
library(ggh4x)


# R version 4.2.2

######## Analysis of groups of neighbourhoods identified by miloR, for more details please see "https://github.com/MarioniLab/miloR" ########
#this script illustrates how we performed our analysis of groups within the endothelial-stromal-cell subset between IBM and CTRL

#load your data

milo_subset <- readRDS(file = "milo_subset.rds")
da_results_subset_IBMvCTRL_adj <- readRDS(file = "da_results_subset_IBMvCTRL_adj.rds")
seed <- readRDS("miloRgroupsSeed.rds")


milo_subset <- buildNhoodGraph(milo_subset)


#calculate groups

set.seed(8)
.Random.seed <- seed

da_results_subset_IBMvCTRL_adj <- groupNhoods(milo_subset, da_results_subset_IBMvCTRL_adj, max.lfc.delta = 3, overlap=2) 


#we identified 16 groups

da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '1'] <- 'Group1'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '2'] <- 'Group2'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '3'] <- 'Group3'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '4'] <- 'Group4'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '5'] <- 'Group5'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '6'] <- 'Group6'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '7'] <- 'Group7'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '8'] <- 'Group8'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '9'] <- 'Group9'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '10'] <- 'Group10'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '11'] <- 'Group11'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '12'] <- 'Group12'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '13'] <- 'Group13'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '14'] <- 'Group14'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '15'] <- 'Group15'
da_results_subset_IBMvCTRL_adj$NhoodGroup[da_results_subset_IBMvCTRL_adj$NhoodGroup == '16'] <- 'Group16'


#create neighborhood group plot
#you can customize the color palette

plotNhoodGroups(milo_subset, da_results_subset_IBMvCTRL_adj, layout="UMAP") +
  scale_fill_manual(values = c("#F512F5", "#F6B0F6", "#FFD500", "#BD1818", "#F9883D", "#A577FF", "#FF8D8D", "#2E008B", "#98BFFA", "#11CDCD", "#438AF5", "#00B7DB", "#1221F7", "#1175CD", "#F75D8C", "#12D0A4")) + 
  force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))


#identify group marker genes
#exclude zero counts genes and find hvgs

keep.rows <- rowSums(logcounts(milo_subset)) != 0
milo_subset <- milo_subset[keep.rows, ]
dec <- modelGeneVar(milo_subset)
hvgs <- getTopHVGs(dec, n=2500)
head(hvgs)  

nhood_markers_subset <- findNhoodGroupMarkers(milo_subset, da_results_subset_IBMvCTRL_adj, subset.row = hvgs, 
                                              aggregate.samples = TRUE, sample_col = "sample_id")

View(nhood_markers_subset)

#use levels to arrange groups for plot creation

da_results_subset_IBMvCTRL_adj$NhoodGroup <- factor(da_results_subset_IBMvCTRL_adj$NhoodGroup, levels = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6", "Group7", "Group8", "Group9", "Group10", "Group11", "Group12", "Group13", "Group14", "Group15", "Group16"))
levels(da_results_subset_IBMvCTRL_adj$NhoodGroup)

#create more plots for visualization

plotDAbeeswarm(da_results_subset_IBMvCTRL_adj, group.by = "NhoodGroup") + 
  force_panelsizes(rows = unit(90, "mm"), cols = unit(120, "mm"))

plotNhoodExpressionGroups(milo_subset, da_results_subset_IBMvCTRL_adj, features= c("insert genes you want to look for"), scale=TRUE,
                          grid.space = "fixed")

