#Summary: Load Vireo results and identify and remove doublets
#Conda: cdh

library(Seurat)
library(doFuture)

registerDoFuture()

set.seed(100)
options(width=100)




donors.dir = "../OutputData/Donors/"
pc.dims = 30

##Filter settings
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

filtered_seurat = readRDS(paste0("../RData/",filter.string,".rds"))


#split object by sample
samples.list = SplitObject(filtered_seurat, split.by = "sample")


#assign Vireo doublets and donors to split samples
sample.res = foreach(sample = names(samples.list)) %do% {
    print(sample)

    sample.barcodes = rownames(samples.list[[sample]]@meta.data)
    
    sample.donors.file = paste0(donors.dir,"/",sample,"/vireo.out/donor_ids.tsv")
    print(sample.donors.file)

    sample.donors = read.table(sample.donors.file, header = T, sep = "\t")
    print(table(sample.donors$best_singlet))
    print(table(sample.donors$donor_id))

    donors = unique(s.donors$best_singlet)
    print(table(s.donors$donor_id))

    sample.doublets = subset(sample.donors, donor_id == "doublet")$cell

    donor.res = foreach(donor = donors) %do% {
        print(donor)

        donor.rows = subset(sample.donors, best_singlet == donor)

        barcode.res = foreach(sample.barcode = sample.barcodes) %do% {

            barcode.parts = strsplit(sample.barcode,"_")[[1]]
            barcode = barcode.parts[length(barcode.parts)]

            if(barcode %in% donor.rows$cell){# is "cell" correct?
                samples.list[[sample]]@meta.data[sample.barcode,"Donor"] = donor
            }

            if(barcode %in% sample.doublets){
                samples.list[[sample]]@meta.data[sample.barcode,"Doublet"] = 1
            } else {
                samples.list[[sample]]@meta.data[sample.barcode,"Doublet"] = 0
            }
        }
    }
}


#combine updated seurat list
filtered_seurat_withdoublets = merge(samples.list[[1]],samples.list[-1])


xtabs(~ sample + Donor, data = filtered_seurat_withdoublets@meta.data, addNA = T)
xtabs(~ sample + Doublet, data = filtered_seurat_withdoublets@meta.data, addNA = T)

#create cell count table, updated with Donors,Doublets
cell.counts = store.counts(so = filtered_seurat_withdoublets, stage = paste0(filter.string,".WithDoublets"), ID.col = "sample")
write.csv(cell.counts, file = paste0("../OutputData/QC/",filter.string,".WithDoublets.csv"), row.names = F)






#create UMAPs of doublet locations

#grab doublet cells for highlighting in UMAP
doublets = rownames(subset(filtered_seurat_withdoublets@meta.data, Doublet == 1)@meta.data)

#need normalized data - will do this in next script as well
#normalize, regress out MT variation, find variable features, scale, run pca
# no established cell-cycle genes for rat, so skipping correction
# return split object
normalized.samples.list = SCT.normalize(filtered_seurat_withdoublets)
filtered_seurat_withdoublets = merge(normalized.samples.list[[1]],normalized.samples.list[-1])
filtered_seurat_withdoublets = RunPCA(filtered_seurat_withdoublets, assay = "SCT")
filtered_seurat_withdoublets = RunUMAP(filtered_seurat_withdoublets, dims = 1:pc.dims,reduction = "pca", return.model = T)
filtered_seurat_withdoublets = RunTSNE(filtered_seurat_withdoublets, check_duplicates = F)

#highlight doublet cells in UMAP
pdf(file=paste0("../Graphs/QC/",filter.string,".Doublets.pdf"))
p = DimPlot(
    filtered_seurat_withdoublets, 
    cells.highlight = doublets,
    split.by = "sample"
)
dev.off()

#remove cells with NAs for Donor assignment - shouldn't be any
dim(filtered_seurat_withdoublets)
filtered_seurat_withdoublets = filtered_seurat_withdoublets[,!is.na(si@meta.data$Donor)]
dim(filtered_seurat_withdoublets)

#remove doublets
filtered_seurat_nodoublets = filtered_seurat_withdoublets[,filtered_seurat_withdoublets@meta.data$Doublet == 0]
dim(filtered_seurat_nodoublets)

table(filtered_seurat_nodoublets$sample)


#redo after removing doublets
#filtered_seurat_nodoublets = RunPCA(filtered_seurat_nodoublets, assay = "SCT")
#filtered_seurat_nodoublets = RunUMAP(filtered_seurat_nodoublets, dims = 1:pc.dims,reduction = "pca", return.model = T)
#filtered_seurat_nodoublets = RunTSNE(filtered_seurat_nodoublets, check_duplicates = F)

#store filtered counts
cell.counts = store.counts(so = filtered_seurat_nodoublets, stage = paste0(filter.string,".NoDoublets"), ID.col = "sample")
write.csv(cell.counts, file = paste0("../OutputData/QC/",filter.string,".NoDoublets.csv"), row.names = F)

#save filtered object
saveRDS(filtered_seurat_nodoublets, file = paste0("../RData/",filter.string,".donors.nodoublets.rds"))


