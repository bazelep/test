#Summary: Integrate 2 E21.5 healthy controls
#Conda: cdh

library(Seurat)
library(ggplot2)

#global ggplot theme
theme_set(theme_bw())

set.seed(100)
options(width=150)

source("utilities.r")


pc.dims = 30

#location to store integration plots
integration.dir = "../Graphs/Integration/"
if(!dir.exists(integration.dir)){dir.create(integration.dir)}


#load data
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
so.list =  readRDS(paste0("../RData/",filter.string,".donors.nodoublets.SCTnormalized_UMIcount,MTratio.rds"))


#grab 2 healthy control E21.5 samples
healthy.samples = c("New_E215_HealthyControl","Original_E215_HealthyControl")

healthys.list = sapply(healthy.samples,function(x){so.list[[x]]})



#first generate a UMAP prior to integration

#combine variable features from each sample
variable.features = unique(
    c(
        healthys.list[[1]]@assays$SCT@var.features,
        healthys.list[[2]]@assays$SCT@var.features
    )
)
length(variable.features)

#merge samples
healthys.preintegration = merge(healthys.list[[1]],healthys.list[-1])


#PCA, UMAP for plot
healthys.preintegration = RunPCA(healthys.preintegration, assay = "SCT", features = variable.features)
healthys.preintegration = RunUMAP(healthys.preintegration, dims = 1:pc.dims,reduction = "pca")
healthys.preintegration = RunTSNE(healthys.preintegration, check_duplicates = F)

pdf(file=paste0(integration.dir,"/","2HealthyControls.PreIntegration.pdf"))
print(DimPlot(healthys.preintegration,reduction = "umap",group.by = "sample"))
dev.off()

#Integration
integ.features = SelectIntegrationFeatures(
    object.list = healthys.list,
    nfeatures = 3000,
    fvf.nfeatures = 3000
) 
length(integ.features)

healthys.list = PrepSCTIntegration(
    object.list = healthys.list,
    anchor.features = integ.features
)

anchors = FindIntegrationAnchors(
    object.list = healthys.list,
    normalization.method = "SCT",
    anchor.features = integ.features,
    k.anchor=5
)

healthys.integrated = IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT"
)

healthys.integrated = RunPCA(healthys.integrated, assay = "integrated")
healthys.integrated = RunUMAP(healthys.integrated, dims = 1:pc.dims,reduction = "pca")
healthys.integrated = RunTSNE(healthys.integrated, check_duplicates = F)

pdf(file=paste0(integration.dir,"/","2HealthyControls.PostIntegration.pdf"))
print(DimPlot(healthys.integrated,reduction = "umap",group.by = "sample"))
dev.off()


saveRDS(healthys.integrated, file = "../RData/2HealthyControls.integrated.rds")
