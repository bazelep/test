#Summary: Add cell type labels to each experimental group object
#Conda: cdh


library(Seurat)
library(read_excel)
library(ggplot2)

#global ggplot theme
theme_set(theme_bw())

set.seed(100)
options(width=100)



clusters.dir = "../OutputData/Clusters"

#load list of seurat objects
so.list = readRDS("../RData/combined.samples.clustered.rds")
statuses = names(so.list)


#load cell type labels table
cell.types.file = "../Data/Cell.Types.xlsx"
cell.types = data.frame(read_excel(cell.types.file, sheet = "Sheet1"))
cell.types = subset(cell.types, !is.na(Cluster))

#create lookup table
#cell.types.list = split(cell.types, f = cell.types$Sample)



status.res = foreach(status = statuses) %do% {

  if(status == "2HC"){
    assay = "integrated"
  } else {
    assay = "SCT"
  }  
  
  so.list[[status]]@meta.data[,"Cluster"] = "NotAssigned"
  so.list[[status]]@meta.data[,"Cell.Type.1"] = "NotAssigned"
  so.list[[status]]@meta.data[,"Cell.Type.2"] = "NotAssigned"
  so.list[[status]]@meta.data[,"Found"] = 0

  cell.type.rows = subset(cell.types, Status == status)
    
  for(cell.type.row.i in c(1:nrow(cell.type.rows))){

    cell.type.row.i.res = cell.type.rows[cell.type.row.i,"Resolution"]
    cell.type.row.i.clu = cell.type.rows[cell.type.row.i,"Cluster"]
    cell.type.row.i.ct1 = cell.type.rows[cell.type.row.i,"Cell.Type.1"]
    cell.type.row.i.ct2 = cell.type.rows[cell.type.row.i,"Cell.Type.2"]

    for(status.row.i in c(1:nrow(so.list[[status]]@meta.data))){
    
      ID.col = paste0(assay,"_snn_res.",cell.type.row.i.res)
    
      status.row.i.cl = so.list[[status]]@meta.data[status.row.i,ID.col]
      
      if(status.row.i.cl == cell.type.row.i.clu){
            
        so.list[[status]]@meta.data[status.row.i,"Cluster"] = paste0(cell.type.row.i.res,".",cell.type.row.i.clu)
        so.list[[status]]@meta.data[status.row.i,"Cell.Type.1"] = cell.type.row.i.ct1
        so.list[[status]]@meta.data[status.row.i,"Cell.Type.2"] = cell.type.row.i.ct2
        so.list[[status]]@meta.data[status.row.i,"Label"] = paste0(cell.type.row.i.res,".",cell.type.row.i.clu,".",cell.type.row.i.ct1)
        so.list[[status]]@meta.data[status.row.i,"Found"] = 1
      }
    }
  }
}



#save table of cell type 1 counts
so = merge(so.list[[1]],so.list[-1])
write.csv(xtabs(~Label+sample,data=so@meta.data), file = paste0(clusters.dir,"/FullLabel.CellCounts.csv"))
write.csv(xtabs(~Cell.Type.1+sample,data=so@meta.data), file = paste0(clusters.dir,"/Cell.Type.1.CellCounts.csv"))
write.csv(xtabs(~Cell.Type.2+sample,data=so@meta.data), file = paste0(clusters.dir,"/Cell.Type.2.CellCounts.csv"))


#plot umap with labels
pdf(file=paste0(clusters.dir,"/2HC,NC,CDH.Labeled.FullLabel.pdf"))
for(s in c("2HC","NC","CDH")){
  print(s)
  Idents(so.list[[s]]) = "Label"
  print(DimPlot(so.list[[s]], label=T, label.size = 2, reduction = "umap",repel=T)+NoLegend()+ggtitle(s))
}
dev.off()

pdf(file=paste0(clusters.dir,"/2HC,NC,CDH.Labeled.Cell.Type.1.pdf"))
for(s in c("2HC","NC","CDH")){
  print(s)
  Idents(so.list[[s]]) = "Cell.Type.1"
  print(DimPlot(so.list[[s]], label=T, label.size = 2, reduction = "umap",repel=T)+NoLegend()+ggtitle(s))
}
dev.off()

pdf(file=paste0(clusters.dir,"/2HC,NC,CDH.Labeled.Cell.Type.2.pdf"))
for(s in c("2HC","NC","CDH")){
  print(s)
  Idents(so.list[[s]]) = "Cell.Type.2"
  print(DimPlot(so.list[[s]], label=T, label.size = 2, reduction = "umap",repel=T)+NoLegend()+ggtitle(s))
}
dev.off()



#save labeled object list
saveRDS(so.list, file = "../RData/2HC,NC,CDH.labeled.list.rds")


