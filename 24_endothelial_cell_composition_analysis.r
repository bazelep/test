#Summary: Compare endothelial cell prportions for pairs of experimental groups.
#Conda: cdh


set.seed(100)
options(width=100)



#output folder
composition.dir = "../OutputData/Composition/"

if(!dir.exists(composition.dir)){
    dir.create(composition.dir)
}


#load labeled data
so.list = readRDS("../RData/2HC,NC,CDH.labeled.list.rds")

#combined seurat objects
so = merge(so.list[[1]],so.list[-1])

#select endothelial cells
Idents(so) = "Cell.Type.2"
so = subset(so, idents = "Endothelial")

#tabulate EC subtypes
xtabs(~Cell.Type.1+Status,data=so@meta.data)

#status pairs
pairs = combn(c("2HC","NC","CDH"),m=2)

#comparing just EC subtypes
idents = c("Cell.Type.1")


for(ident in idents){
    results = data.frame()

    print(ident)
    
    Idents(so) = ident
    
    #get metadata
    metadata.cols = c("sample","Donor",ident,"Status")
    metadata = so@meta.data[,metadata.cols]

    #tabulate EC subtypes by Status
    # - each row is an EC subtype
    # - each column is a status
    xt = xtabs(as.formula(paste0("~",ident,"+Status")),data=metadata)
    
    #iterate over pairs of statuses
    for(pair.i in c(1:ncol(pairs))){

        pair = pairs[,pair.i]
        print(p)
        
        xt.pair = xt[,pair]
        print(xt.pair)
        
        total.cell.counts = colSums(xt.pair)
        
        #iterate over cell types
        for(cell.type in rownames(xt.pair)){
            print(cell.type)

            #grab counts for cell type
            cell.type.cell.counts = xt.pair[cell.type,]
            print(cell.type.cell.counts)

            #if >1 cell for each type and status
            if(all(cell.type.cell.counts > 0)){

                #compare proportions of cell type for each status in pair
                fit = prop.test(x=cell.type.cell.counts,n=total.cell.counts,alternative="two.sided",correct=F)

                #create result table
                res = data.frame(
                    "Cell.Type.Grouping"=ident,
                    "Cell.Type"=cell.type,
                    "Group1"=pair[1],
                    "Group2"=pair[2],
                    "Group1.CellType.CellCounts"=cell.type.cell.counts[1],
                    "Group2.CellType.CellCounts"=cell.type.cell.counts[2],
                    "Group1.Total.CellCounts"=total.cell.counts[1],
                    "Group2.Total.CellCounts"=total.cell.counts[2],
                    "Group1.Proportion"=fit$estimate[1],
                    "Group2.Proportion"=fit$estimate[2],
                    "Pearson.Chisquared.Statistic"=fit$statistic,
                    "P.Value"=fit$p.value
                    )
                results = rbind(results,res)
            }
        }
    }
    results$BH.P.Value = p.adjust(results$P.Value,method="BH")
    results$Bonferroni.P.Value = p.adjust(results$P.Value,method="bonferroni")
    write.csv(results, file = paste0(composition.dir,"/2HC,NC,CDH.labeled.list.ECsSubtypes.",ident,".CompositionAnalysis.csv"),row.names=F)
}





