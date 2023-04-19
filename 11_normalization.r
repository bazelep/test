#Summary: Normalize counts and correct for nuisance variation
#Conda: cdh

library(Seurat)


set.seed(100)
options(width=150)

source("utilities.r")


#Steps:
#1. Normalize for UMI counts per cell (sequencing depth)
#2. Find highly variable genes
#3. Scale data for PCA
#4. PCA



#set filter string for filtered seurat object
f.genes = 500
f.umi = 500
f.mt = 0.2
f.complexity = 0.8
f.rb = 0.2

filter.string = paste0(
    "Filtered.G",f.genes,
    ".U",f.umi,
    ".M",f.mt,
    ".R",f.rb,
    ".C",f.complexity
)


#load data
filtered_seurat = readRDS(file = paste0("../RData/",filter.string,".donors.nodoublets.rds"))

filtered_seurat

#normalize, regress out MT variation, find variable features, scale, run pca
# no established cell-cycle genes for rat, so skipping correction
# return split object, which is needed for integration
so.list = SCT.normalize(filtered_seurat)

#save split object
saveRDS(
    so.list, 
    file = paste0("../RData/",filter.string,".donors.nodoublets.SCTnormalized_UMIcount,MTratio.rds")
)






