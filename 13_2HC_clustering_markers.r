#Summary: Generate clusters and markers for 2 integrated E21.5 healthy controls
#Conda: cdh

library(AnnotationHub)
library(Seurat)
library(rlist)
library(ggplot2)

#global ggplot theme
theme_set(theme_bw())

set.seed(100)
options(width=100)

source("utilities.r")



#get gene descriptions to add to marker genes
species = "Rattus_norvegicus"
snapshot.date = "2021-10-20"

annotations = get.annotations(species = species, snapshot.date = snapshot.date)
head(annotations)


#load integrated data
healthys.integrated = readRDS("../RData/2HealthyControls.integrated.rds")
head(healthys.integrated@meta.data)


#generate SNN graph and clusters for a range of resolutions
resolutions = seq(0,1,by=0.05)
resolutions

healthys.integrated = make.clusters(
    so = healthys.integrated, 
    resolutions = resolutions, 
    sample.set = "E215.2HealthyControls",
    sample.col = "sample"
)

#generate marker genes for each resolution's clusters
#clusters are skipped if not present in all samples/groups
for(res in resolutions){
  print(paste0("...looking at resolution ", res))
  
  get.markers(
    so = healthys.integrated, 
    ID.col = paste0("integrated_snn_res.",res),
    annotations = annotations,
    sample.set = "E215.2HealthyControls",
    group.samples.by = "sample"
  )
}


saveRDS(healthys.integrated, file = "../RData/2HealthyControls.integrated.clustered.rds")
