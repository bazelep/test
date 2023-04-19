#Summary: Split up receptor-ligand comparisons with liana into separate SingleCellExperiment objects to faciliate parallel analysis.
#Conda: cdh


library(Seurat)
library(readxl)
library(doFuture)

set.seed(100)
options(width=100)



#output folder
ligand.receptor.dir = "../OutputData/LigandReceptor/"
comparisons.dir = paste0(ligand.receptor.dir, "/Comparisons")

if(!dir.exists(comparisons.dir)){
    dir.create(comparisons.dir)
}

#load rat to human ortholog map
features.orthologs = readRDS(paste0(ligand.receptor.dir,"/rat_to_human_orthologs.rds"))

#exclude genes of length 0
features.orthologs = subset(features.orthologs, nchar(Human.gene.name) != 0)
features.orthologs = subset(features.orthologs, nchar(Rat.gene.name) != 0)





#load labeled data
so.list = readRDS("../RData/2HC,NC,CDH.labeled.list.rds")


#perivascular fibroblasts has 2 clusters in 2HC which should be combined
xtabs(~Cell.Type.1+sample,data=so.list[["2HC"]]@meta.data)

so.list[["2HC"]]@meta.data[so.list[["2HC"]]$Cell.Type.1 == "Perivascular Fibroblasts Additional","Cell.Type.1"] = "Perivascular Fibroblasts"

xtabs(~Cell.Type.1+sample,data=so.list[["2HC"]]@meta.data)




#load ligand-receptor cell type pairs to compare
comparisons = data.frame(read_excel("../Data/Ligand.Receptor.Comparisons.xlsx", sheet = "Sheet1"))




#note that each status was created separately so may have different genes
lapply(so.list,function(x){dim(x@assays$RNA@counts)})


#find genes found in all 3 statuses

#find all genes
genes = c()
genes.res = foreach(status = names(so.list)) %do% {
  print(status)
  so = so.list[[status]]
  genes = c(genes, row.names(so@assays$RNA@counts))
}
unique.genes = unique(genes)

#filter for present in all 3 statuses
# - create vector of T/F
unique.res = foreach(gene = unique.genes) %do% {
  string = paste0("'^",gene,"$'")
  length(grep(eval(parse(text=string)),genes,value=T)) == 3
}
names(unique.res) = unique.genes
unique.res = unlist(unique.res)

#grab TRUEs
unique.genes.keep = unique.genes[unique.res]
length(unique.genes.keep)
length(unique.genes)

#make sure all have length > 0
range(nchar(unique.genes.keep))


#loop over statuses and grab selected genes
so.list.orthologs = foreach(status = names(so.list)) %do% {
  
    so = so.list[[status]]
    counts = as.matrix(so@assays$RNA@counts)
  
    #keep selected genes
    counts = counts[unique.genes.keep,]

    #keep genes found in ortholog map
    counts = counts[rownames(counts) %in% features.orthologs$Rat.gene.name,]
    
    #extract ortholog map entries found in counts           
    counts.orthologs = features.orthologs[features.orthologs$Rat.gene.name %in% rownames(counts),]
    
    #remove duplicate entries (e.g. ensembl ID maps to >1 gene symbol)
    counts.orthologs = counts.orthologs[!duplicated(counts.orthologs$Rat.gene.name),]
    
    #order counts by gene symbol so can replace rat gene symbol with human gene symbol
    counts = counts[counts.orthologs$Rat.gene.name,]

    #update count gene names (rows) with human ortholog
    rownames(counts) = counts.orthologs$Human.gene.name

    #create new seurat object with modified counts - metadata doesn't change
    new.so = CreateSeuratObject(counts = counts, meta.data = so@meta.data)

    return(new.so)
}
names(so.list.orthologs) = names(so.list)

#check for new gene counts to make sure are same for all 3 statuses
lapply(so.list.orthologs,function(x){dim(x@assays$RNA@counts)})



#vector to track comparisons with missing ligand/receptor cell types
missing.comps = c()

#iterate over each comparison cell type pair
comps.res = foreach(comp.i = c(1:nrow(comparisons))) %do% {

    #list to store ligand and receptor cells
    cells = list()

    #grab a cell type pair to be compared 
    comp.row = comps[comp.i,]
    print(comp.row)

    #iterate over each status
    samples.res = foreach(status = names(so.list.orthologs)) %do% {
    
        print(status)

        #if cell type missing, will flag
        missing = F

        #grab seurat object
        so = so.list.orthologs[[s]]

        #apply Ca4+ expression filter if present
        if(comp.row$Cell.Subset == "Ca4plus"){
            so = subset(so, Ca4plus == T)
        }

        #extract ligand cells and then receptor cells by alternately setting identities of so

        #assign ligand cell type label level ('Cell.Type.1', or 'Cell.Type.2')
        Idents(so) = comp.row$Ligand.Group.Variable
        print(table(Idents(so)))
        
        #if ligand cell group present, subset
        if (comp.row$Ligand.Group %in% Idents(so)){
            cells[["ligand"]] = subset(so, idents = comp.row$Ligand.Group)
        } else {
            missing = T
        }
        
        #assign receptor cell type label level ('Cell.Type.1', or 'Cell.Type.2')
        Idents(so) = comp.row$Receptor.Group.Variable
        print(table(Idents(so)))
        
        #if receptor cell group present, subset
        if (comp.row$Receptor.Group %in% Idents(so)){
            cells[["receptor"]] = subset(so, idents = comp.row$Receptor.Group)
        } else {
            missing = T
        }

        #if neither ligand nor receptor cell type missing, proceed
        if(!missing){        

            #create list of SingleCellExperiment objects for liana
            sce.list = sapply(names(cells), function(cell.group){

                #add metadata label for inclusion in SCE object
                cells[[cell.group]]$Ligand.Receptor.Group = cell.group

                #grab cell counts
                cell.group.counts = as.matrix(cells[[cell.group]]@assays$RNA@counts)
                
                #grab metadata
                cell.group.coldata = cells[[cell.group]]@meta.data[,c("cells","Ligand.Receptor.Group")]
                
                #create SCE object
                cell.group.sce = SingleCellExperiment(
                  assays = list(counts=cell.group.counts),
                  colData = cell.group.coldata
                )
                
                return(cell.group.sce)
            })

            #combined ligand and receptor SCE objects
            sce = cbind(sce.list[[1]],sce.list[[2]])

            #check column (cell) numbers
            print(ncol(sce))

            #save SCE object, named by comparison row number and status
            save(sce, file = paste0(comparisons.dir,"/",comp.row$Comparison,".",status,".sce.RData"))      

        } else { #one/both cell types was missing
          missing.comps = c(missing.comps, comp.row$Comparison)
          print(paste0("..missing identity for comp ",comp.row$Comparison,"..."))
        }
    }
}

#any missing?
comparisons[missing.comps,]


