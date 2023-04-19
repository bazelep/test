#Summary: Load counts, create seurat objects, add QC metrics
#Conda: cdh

library(Seurat)
library(doFuture)


registerDoFuture()
plan("sequential")

set.seed(100)




options(width=150)




counts.base.dir = "../OutputData/Counts/"



#load cellranger counts and create initial seurat objects
counts.dirs = list.files(path = counts.base.dir, full.names=T)
counts.dirs

counts.dirs = grep("mro",grep("log",counts.dirs,invert=T,value=T),invert=T,value=T)


samples.list = foreach(counts.dir = counts.dirs) %do% {
    print(counts.dir)
    
    matrix.dir = paste0(counts.dir,"/outs/filtered_feature_bc_matrix/")
    print(matrix.dir)

    s.mat = Read10X(data.dir = matrix.dir)

    # Initialize the Seurat object 
    so = CreateSeuratObject(
      counts = s.mat,
      project = basename(counts.dir), 
      min.cells = 1, 
      min.features = 1
    )
                               
    print(dim(so))
    return(so)
}
names(samples.list) = basename(counts.dirs)
samples.list



#Not including E14.5 or E18.5 anymore, or original E21.5 Nitrofen
remove = c(
  grep("E14_5",names(samples.list)),
  grep("E18.5",names(samples.list)),
  grep("E21_5_NitrofenFetalRatLung",names(samples.list))
)

samples.list = samples.list[-remove]
names(samples.list)


#update sample names for newer set of samples
sample.map = list()
sample.map[["Sample_FC_JA_2998_1"]] = "New_E215_HealthyControl"
sample.map[["Sample_FC_JR_3032_1"]] = "E215_NitrofenControl"
sample.map[["Sample_FC_JR_3032_2"]] = "E215_NitrofenCDH"
sample.map[["E21_5_LeftLungControl"]] = "Original_E215_HealthyControl"



#check original sample counts
t(t(sapply(1:length(samples.list), function(x){table(samples.list[[x]]$orig.ident)})))


#add 'sample' metadata column containing original or updated sample name
map.res = foreach(sample = names(samples.list)) %do% {
  print(sample)

  samples.list[[sample]]$sample = samples.list[[sample]]$orig.ident
  print(xtabs(~orig.ident+sample,data=samples.list[[sample]]@meta.data))

  original.name = as.character(samples.list[[sample]]$orig.ident[1])
  print(original.name)

  if(original.name %in% names(sample.map)){
    print(sample.map[[original.name]])
    samples.list[[sample]]$sample = sample.map[[original.name]]
  }
}

#check updated names
t(t(sapply(1:length(samples.list), function(x){table(samples.list[[x]]$sample)})))


saveRDS(samples.list, file = "../RData/samples.list.rds")


#merge unfiltered samples
merged_seurat = merge(x=samples.list[[1]],y=samples.list[-1],add.cell.id = names(samples.list))

head(merged_seurat@meta.data)

table(merged_seurat@meta.data$orig.ident)
xtabs(~orig.ident+sample,data=merged_seurat@meta.data)


##Add QC Metrics

#Complexity
merged_seurat$log10GenesPerUMI = log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

#Mitochondrial genes proportion
merged_seurat$mitoRatio = PercentageFeatureSet(object = merged_seurat, pattern = "^Mt-") / 100

#Ribosomal genes proportion
cts = GetAssayData(merged_seurat, assay = "RNA", slot = "counts")
rb.genes = grep("^RP[SL]",rownames(cts),value=T)
merged_seurat$riboRatio = colSums(cts[rb.genes,])/colSums(cts)
rm(cts)


#add for convenience
merged_seurat$barcode = rownames(merged_seurat@meta.data)


saveRDS(merged_seurat, file = "../RData/merged_seurat.rds")

