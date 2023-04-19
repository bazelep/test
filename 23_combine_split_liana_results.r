#Summary: Combine results from split up receptor-ligand comparisons with liana.
#Conda: cdh

library(readxl)
library(doFuture)

set.seed(100)
options(width=100)




#output folder
ligand.receptor.dir = "../OutputData/LigandReceptor/"
comparisons.dir = paste0(ligand.receptor.dir, "/Comparisons")


#load cell type comparison pairs
comparisons = data.frame(read_excel("../Data/Receptor.Ligand.Comparisons2.xlsx", sheet = "Sheet1"))


#get list of liana result files
result.files = list.files(path=comparisons.dir,pattern="results\\.RData",full.names=T)
result.files



registerDoFuture()
options(future.globals.maxSize = 8000 * 1024^2)
plan("sequential")


#iterate over result files and gather consensus/aggregate rank results
files.res = foreach(file = result.files) %do% {

    comp = strsplit(basename(file),"\\.")[[1]][1]
    status = strsplit(basename(file),"\\.")[[1]][2]

    liana.results = readRDS(file)

    agg = data.frame(liana.results[["aggregate"]])
    agg$Comparison = comp
    agg$Status = status

    return(agg)
}
agg.results = list.rbind(files.res)
dim(agg.results)
head(agg.results)

#liana compares both cell types as both receptor and ligand, so there are duplicate results, so remove these
keep = subset(agg.results, (source != target) & (source == "ligand") & (ligand.complex != receptor.complex) )

#combine results with comparisons table
comps.keep = merge(comps, keep, by = "Comparison")

#order results by aggregate.rank column 
comps.keep = comps.keep[order(comps.keep$aggregate_rank),]

#how many significant?
sum(comps.keep$aggregate_rank < 0.05)

#save results
#save(comps.keep, file = paste0(comparisons.dir,"/liana.consensus.results.pergene.RData"))
write.csv(comps.keep, file = paste0(comparisons.dir,"/liana.consensus.results.pergene.csv"))



#create pairwise matrices of counts of significant (aggregate rank < 0.05) genes
# - each row is a ligand cell type
# - each column is a receptor cell type
# - separate matrices for each status

#matrix cell types to include
cell.types = c(
"Endothelial",
"mv ECs",
"mv Proliferative ECs",
"High Hg ECs",
"mv ECs high in Ca4",
"Epithelial",
"AT1",
"AT2",
"BPs",
"AllEC.Ca4plus",
"mvEC.Ca4plus",
"Matrix Fibroblasts",
"Perivascular Fibroblasts",
"Vascular SMCs",
"Mature Pericytes",
"SMC-like Pericytes",
"Neutrophils"
)

statuses = c("2HC","NC","CDH")

#get significant genes
comps.keep.sig = subset(comps.keep, aggregate_rank < 0.05)
dim(comps.keep)
dim(comps.keep.sig)


#create emptry list for each status and cell type, which will be used to construct matrices
mat.list = list()
for(status in statuses){
    for(cell.type.i in cell.types){
        for(cell.type.j in cell.types){
            mat.list[[status]][[cell.type.i]][[cell.type.j]] = c()
        }
    }
}

#iterate over significant genes and append to list
for(row in c(1:nrow(comps.keep.sig))){

    status = comps.keep.sig[row,"Status"]

    ligand.cell.type = comps.keep.sig[row,"Ligand.Group"]
    
    receptor.cell.type = comps.keep.sig[row,"Receptor.Group"]

    #create string of ligand and receptor genes
    pair = paste0(
        comps.keep.sig[row,"ligand.complex"],
        ";",
        comps.keep.sig[row,"receptor.complex"]
    )

    #append gene pair
    mat.list[[status]][[ligand.cell.type]][[receptor.cell.type]] = c(
        mat.list[[status]][[ligand.cell.type]][[receptor.cell.type]], 
        pair
    )
}



#iterate over statuses and create matrices
for(status in c("2HC","NC","CDH")){
  
    #create empty matrix
    mat = matrix(0, nrow = length(cell.types), ncol = length(cell.types))
    rownames(mat) = cell.types
    colnames(mat) = cell.types

    #iterate over ligand cell types (rows of matrix)
    for(cell.type.i in c(1:length(cell.types))){
        print(cell.type.i)

        #iterate over receptor cell types (columns of matrix)
        for(cell.type.j in c(1:length(cell.types))){
            print(cell.type.j)

            #grab cell type labels
            cell.type.i.label = cell.types[cell.type.i]
            cell.type.j.label = cell.types[cell.type.j]
            print(cell.type.i.label)
            print(cell.type.j.label)

            #initialize gene count with NA in case no entries present
            significant.gene.count = NA

            #if ligand cell type present for status
            if(cell.type.i.label %in% names(mat.list[[status]])){
      
                #if receptor cell type present for status
                if(cell.type.j.label %in% names(mat.list[[s]][[cell.type.i.label]])){
      
                    #get count of significant genes for that ligand-receptor pair
                    significant.gene.count = length(unique(mat.list[[status]][[cell.type.i.label]][[cell.type.j.label]]))
                }
            } 
      
            print(significant.gene.count)

            #add gene count to matrix
            mat[cell.type.i,cell.type.j] = significant.gene.count
        }
    }
    print(mat)
    
    #save matrix
    write.csv(mat, file = paste0(comparisons.dir,"/Liana.AggregateRank.LessThan0.05.",status,".Matrix.csv"))
}



#gather non-aggregate/consensus results from liana (individual method results)
# - default methods were included
# - only cellphonedb calculates p-values, so only filtered these

plan("multisession", workers = 20)


#iterate over files
files.res = foreach(file = result.files) %dopar% {

    comp = strsplit(basename(file),"\\.")[[1]][1]
    status = strsplit(basename(file),"\\.")[[1]][2]

    #load liana results rds
    liana.results = readRDS(file)

    sep = data.frame(liana.results[["separate"]])

    #names of methods used by liana
    methods = unique(sapply(strsplit(names(sep),"\\."),"[",1L))

    #iterate over methods
    method.res = foreach(method = methods) %dopar% {
      #print(method)

        #grab method-specific columns
        method.cols = sep[,grep(method,names(sep),value=T)]

        #remove method names from columns
        names(method.cols) = gsub(paste0(method,"\\."),"",names(method.cols))

        #add annotation
        method.cols$Method = method
        method.cols$Comparison = comp
        method.cols$Status = status

        #remove duplicates as for aggregat results
        method.cols.keep = subset(method.cols, (source != target) & (source == "ligand") & (ligand.complex != receptor.complex) )
      
        #combine with comparisons table
        comp.method.cols.keep = merge(comps, method.cols.keep, by = "Comparison")
      
        #order columns based on interaction metric
        if(method == 'natmi'){
            comp.method.cols.keep = comp.method.cols.keep[order(comp.method.cols.keep$edge_specificity,decreasing=T),]
            comp.method.cols.keep.sig = comp.method.cols.keep
        } else if(method == 'connectome'){
            comp.method.cols.keep = comp.method.cols.keep[order(comp.method.cols.keep$weight_sc,decreasing=T),]
            comp.method.cols.keep.sig = comp.method.cols.keep
        } else if(method == 'logfc'){
            comp.method.cols.keep = comp.method.cols.keep[order(comp.method.cols.keep$logfc_comb,decreasing=T),]
            comp.method.cols.keep.sig = comp.method.cols.keep
        } else if(method == 'sca'){
            comp.method.cols.keep = comp.method.cols.keep[order(comp.method.cols.keep$LRscore,decreasing=T),]
            comp.method.cols.keep.sig = comp.method.cols.keep
        } else { #cellphonedb
            comp.method.cols.keep = comp.method.cols.keep[order(comp.method.cols.keep$pvalue),]

            #filter by p-value
            comp.method.cols.keep.sig = subset(comp.method.cols.keep, pvalue < 0.05)
        }
  
        return(comp.method.cols.keep.sig)
    }
    names(method.res) = methods    
    return(method.res)
}
names(files.res) = basename(result.files)


#save results for each method separately, since have different columns
methods.list = list()

#need to invert methods and comparisons so can aggregate by method
list.res = foreach(comp = names(files.res)) %do% {
  
    #iterate over method
    method.res = foreach(method = names(files.res[[comp]])) %do% {

        #create inverted list
        methods.list[[method]][[comp]] = files.res[[comp]][[method]]
    }
}

#save results for each method
methods.res = foreach(method = names(methods.list)) %do% {

  res = list.rbind(methods.list[[method]])
  write.csv(res, file = paste0(comparisons.dir,"/liana.",method,".results.pergene.csv"),row.names=F)
}


