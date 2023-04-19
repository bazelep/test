#Summary: Generate per-gene expression percentages and z-scores for dot plots
#Conda: cdh


library(Seurat)

set.seed(100)
options(width=100)



#output folder
markers.dir = "../OutputData/MarkerGenes/"




#load labeled data
so.list = readRDS("../RData/2HC,NC,CDH.labeled.list.rds")


#perivascular fibroblasts has 2 clusters in 2HC which should be combined
xtabs(~Cell.Type.1+sample,data=so.list[["2HC"]]@meta.data)

so.list[["2HC"]]@meta.data[so.list[["2HC"]]$Cell.Type.1 == "Perivascular Fibroblasts Additional","Cell.Type.1"] = "Perivascular Fibroblasts"

xtabs(~Cell.Type.1+sample,data=so.list[["2HC"]]@meta.data)



#calculate at both granular and global cell type assignments
idents = c("Cell.Type.1","Cell.Type.2")


#calculate percentage of cells expression a gene, for each gene

#calculate mean library-size-normalized expression across cells for a gene, for each gene
# - use NormalizeData normalized data - need to take exponential and add 1 before averaging


for(ident in idents){
    print(ident)

    #output folder
    res.dir = paste0(markers.dir,"/",ident,"/SMCs")
    
    if(!dir.exists(res.dir)){
        dir.create(res.dir,recursive=T)
    }
    
    ident.results = data.frame()
    
    #iterate over experimental group
    for(status in c("2HC","NC","CDH")){
        print(status)

        Idents(so.list[[status]]) = ident
        
        #iterate over clusters
        for(cl in clusters){
            print(cl)

            #subset to cluster
            so.cl = subset(so.list[[s]], idents = cl)      

            #grab raw counts
            cl.cts = GetAssayData(so.cl, assay = "RNA", slot = "counts")

            #grab NormalizeData normalized counts
            cl.data = GetAssayData(so.cl, assay = "RNA", slot = "data")

            expr.pct = function(g){
                return(sum(cl.cts[g,] > 0)/ ncol(cl.cts))
            }
            expr.avg = function(g){
                return(mean(expm1(cl.data[g,])))
            }

            #calculate expression percentages for each gene
            pcts = data.frame(sapply(rownames(cl.cts), expr.pct))
            
            #calculate expression averages for each gene
            avgs = data.frame(sapply(rownames(cl.cts), expr.avg))
            
            #merge based on genes
            pcts.avgs = merge(pcts,avgs,by="row.names")
            
            #rename columns
            names(pcts.avgs) = c("Gene","pct.exp","avg.exp")

            #scale across genes
            max = max(pcts.avgs$avg.exp)
            min = min(pcts.avgs$avg.exp)
            
            #calculate z-score for each gene
            pcts.avgs$avg.exp.scaled.zscore = (pcts.avgs$avg.exp - mean(pcts.avgs$avg.exp)) / sd(pcts.avgs$avg.exp)

            #set high and low numbers to 2.5 - the rationale is Seurat's DotPlot function
            pcts.avgs$avg.exp.scaled.zscore.cut = ifelse(pcts.avgs$avg.exp.scaled.zscore > 2.5, 2.5, pcts.avgs$avg.exp.scaled.zscore)
            pcts.avgs$avg.exp.scaled.zscore.cut = ifelse(pcts.avgs$avg.exp.scaled.zscore.cut < -2.5, -2.5, pcts.avgs$avg.exp.scaled.zscore.cut)
            
            #scale from 0 to 1
            pcts.avgs$avg.exp.scaled.range01 = (pcts.avgs$avg.exp - min) / (max - min)
            
            #add meta-data
            pcts.avgs$Cluster = cl
            pcts.avgs$Cell.Type.Group = ident
            pcts.avgs$Status = status
            

            #change column order of output

            #cell column index
            ct.pos = grep("Cell",colnames(pcts.avgs))
            
            #status column index
            s.pos = grep("Status",colnames(pcts.avgs))
            
            #cluster column index
            clu.pos = grep("Cluster",colnames(pcts.avgs))
            
            #gene column index
            gene.pos = grep("Gene",colnames(pcts.avgs))
            
            #combine ordering of above columns
            first = c(ct.pos,clu.pos,s.pos,gene.pos)
            
            #full set of indices
            cols = c(1:ncol(pcts.avgs))
            
            #rest of the columns
            other.cols = cols[!cols %in% first]
            
            #new ordering
            new.cols = c(first, other.cols)
            
            #reorder
            pcts.avgs = pcts.avgs[,new.cols]  
            
            #combine per-cluster results
            ident.results = rbind(ident.results, pcts.avgs)
        }
    }
    write.csv(ident.results, file = paste0(res.dir,"/",ident,".ExpressionPercentages.csv"), row.names = F)
}





