
#Summary: Apply QC filtering
#Conda: cdh

library(Seurat)
library(ggplot2)

#global ggplot theme
theme_set(theme_bw())

set.seed(100)
options(width=150)

source("utilities.r")



#folders for QC results
qc.counts.dir = "../OutputData/QC/"
qc.plots.dir = "../Graphs/QC/"

if(!dir.exists(qc.counts.dir)){dir.create(qc.counts.dir)}
if(!dir.exists(qc.plots.dir)){dir.create(qc.plots.dir)}



#load merged data
merged_seurat = readRDS(file="../RData/merged_seurat.rds")
merged_seurat

Idents(merged_seurat) = "sample"
table(Idents(merged_seurat))


#Create initial QC plots/cell count table using default filter criteria

##Filter settings
f.genes = 500
f.umi = 500
f.mt = 0.2
f.complexity = 0.8
f.rb = 0.2


filter.string = paste0(
  "Unfiltered.G",f.genes,
  ".U",f.umi,
  ".M",f.mt,
  ".R",f.rb,
  ".C",f.complexity
)

#store initial counts to csv
cell.counts = store.counts(merged_seurat, stage = filter.string)
write.csv(cell.counts, file = paste0(qc.counts.dir,filter.string,".csv"))


qc.plots(
  merged_seurat,
  stage = filter.string,
  f.genes = list("x.min" = NULL,"x.max" = NULL, "x.intercept" = f.genes),
  f.umi = list("x.min" = NULL, "x.max" = NULL, "x.intercept" = f.umi),
  f.mt = f.mt,
  f.rb = f.rb,
  f.complexity = f.complexity,
  qc.dir = qc.plots.dir,
  ID.col = "sample"
)




##Apply Filtering


#Cell-level filtering

# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat = subset(x = merged_seurat, 
                          subset = (nCount_RNA >= f.umi) & 
                            (nFeature_RNA >= f.genes) & 
                            (log10GenesPerUMI > f.complexity) & 
                            (mitoRatio < f.mt))


#Gene-level filtering
# remove genes if expressed in <10 cells

counts = GetAssayData(filtered_seurat, assay = "RNA", slot = "counts")

#create matrix of expressed (>0) 
nonzero = counts > 0
sum(nonzero)

#create vector of genes expressed in >= 10 cells
keep_genes = rowSums(nonzero) >= 10
length(keep_genes)
sum(keep_genes)

#filter counts
filtered_counts = counts[keep_genes, ]

#reassign to filtered Seurat object
filtered_seurat = CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
Idents(filtered_seurat) = "sample"



table(Idents(merged_seurat))
table(Idents(filtered_seurat))

#Create updated table/plots post-filtering
filter.string = paste0(
    "Filtered.G",f.genes,
    ".U",f.umi,
    ".M",f.mt,
    ".R",f.rb,
    ".C",f.complexity
)

#store filtered cell counts
cell.counts = store.counts(filtered_seurat, stage = filter.string)
write.csv(cell.counts, file = paste0(qc.counts.dir,filter.string,".csv"))

#create filtered plots
qc.plots(
  filtered_seurat,
  stage = filter.string,
  f.genes = list("x.min" = NULL,"x.max" = NULL, "x.intercept" = f.genes),
  f.umi = list("x.min" = NULL, "x.max" = NULL, "x.intercept" = f.umi),
  f.mt = f.mt,
  f.rb = f.rb,
  f.complexity = f.complexity,
  qc.dir = qc.plots.dir,
  ID.col = "sample"
)


#save filtered object
saveRDS(filtered_seurat, file = paste0("../RData/",filter.string,".rds"))

