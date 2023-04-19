#Summary: Generate clusters for NC and CDH
#Conda: cdh


library(Seurat)
library(ggplot2)

#global ggplot theme
theme_set(theme_bw())

set.seed(100)
options(width=100)

source("utilities.r")


resolutions = seq(0,1,by=0.05)
resolutions
pc.dims = 30


#get gene descriptions to add to marker genes
species = "Rattus_norvegicus"
snapshot.date = "2021-10-20"

annotations = get.annotations(species = species, snapshot.date = snapshot.date)
head(annotations)



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

nitrofen.samples = c("E215_NitrofenControl","E215_NitrofenCDH")


sample.res = foreach(sample = nitrofen.samples) %do% {
    print(sample)
  
    so.list[[sample]] = RunUMAP(so.list[[sample]], dims = 1:pc.dims,reduction = "pca", return.model = T)
  
    so.list[[sample]] = make.clusters(
        so = so.list[[sample]], 
        resolutions = resolutions, 
        sample.set = sample, 
        sample.col = "sample",
        assay = "SCT"
    )
  
    resolution.res = foreach(res = resolutions) %do% {
        print(paste0("...looking at resolution ",res))
        
        get.markers(
            so = so.list[[sample]], 
            resolution = res,
            ID.col = paste0("SCT_snn_res.",res),
            annotations = annotations,
            sample.set = sample,
            group.samples.by = "sample"
        )
    }

    return(so.list[[sample]])
}
names(sample.res) = nitrofen.samples

#save objects
saveRDS(sample.res, file = "../RData/2NitrofenSamples.clustered.rds")







