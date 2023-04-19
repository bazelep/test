#Summary: Generate pseudobulk counts and perform differential expression testing for experimental group comparisons
#Conda: cdh


library(Seurat)
library(Matrix.utils)
library(doFuture)
library(rlist)

set.seed(100)
options(width=100)




registerDoFuture()
options(future.globals.maxSize = 8000 * 1024^2)
plan("multisession", workers = 30)




#get gene annotation
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <- useEnsembl(biomart = "ensembl", version = 108)
ensembl.rat <- useDataset(dataset = "rnorvegicus_gene_ensembl", mart = ensembl)
rat.anno = getBM(
      attributes= c("external_gene_name","description"),
      mart=ensembl.rat)
rat.anno = rat.anno[order(rat.anno$description,decreasing=T),]

#remove duplicates by gene symbol
dim(rat.anno)
rat.anno = rat.anno[!duplicated(rat.anno$external_gene_name),]
dim(rat.anno)



#experimental group contrasts
coefs = c(
  "CDH - NC",
  "CDH - HC",
  "NC - HC"
)




#limma wrapper
de.limma = function(
    counts = NULL, #pseudobulk count matrix
    samples = NULL, #sample phenotype metadata
    res.dir = "../OutputData/DGE", #location for output
    cluster = NULL #cluster label for output
){

    res.file.base = paste0(res.dir,"/Limma/",cluster)

    if(!dir.exists(res.file.base)){
        dir.create(res.file.base, recursive = T)
    }
    
    #create dgelist
    x = DGEList(counts = counts, samples = samples)

    #TMM normalization
    y = calcNormFactors(x)
    print("...initial normalization factors:")
    print(y$samples$norm.factors) 

    #create design matrix
    design = model.matrix(~0+Status,data=y$samples)
    colnames(design) = gsub("Status","",colnames(design))
    print(design)

    #filter for lowly express genes
    keep.exprs = filterByExpr(y, design = design)

    #re-normalize after removing lowly-expressed genes
    y = calcNormFactors(y[keep.exprs,,keep.lib.sizes=FALSE])
    print("...normalization factors after filtering:")
    print(y$samples$norm.factors)

    #creae mean-difference plots for each sample
    res.file.meandiff = paste0(res.file.base,"/meandiff.pdf")
    pdf(file = res.file.meandiff)
    par(mfrow=c(2,3))
    for (i in seq_len(ncol(y))) {
        plotMD(y, column=i)
    }
    dev.off()

    #calculate log-CPM
    y.lcpm = edgeR::cpm(y, log = T)
    lognorm = y.lcpm
    colnames(lognorm) = paste0(colnames(lognorm), ".TMMLog2CPM")
    write.csv(lognorm, file = paste0(res.file.base,"/","Cluster",cluster,".TMM.Log2CPM.csv"))

    #create MDS plots
    pdf(file = paste0(res.file.base, "/MDS.Status.pdf"))
    plotMDS(y.lcpm, col=y$samples$Color)
    plotMDS(y.lcpm, col=y$samples$Color, dim.plot = c(3,4))
    dev.off()    
    
    #calculate voom precision weights
    v <- voom(y, design = design, plot = F)

    #fit gene models
    fit <- lmFit(v, design = design)
    
    #construct list of contrasts
    myargs = list()
    for(coef in coefs){
        myargs = list.append(myargs, coef)
    }
    myargs = list.append(myargs,levels=design)

    #calculate all contrasts
    contrasts = do.call(makeContrasts, myargs)
    con.fit <- contrasts.fit(fit, contrasts)    

    #ebayes step
    efit <- eBayes(con.fit)  

    #combine results for each contrast    
    res.df = data.frame()
    for(coef in coefs){
        print(coef)

        tt = topTable(efit, coef = coef, sort.by = "P", n = Inf)
        
        tt$Gene = rownames(tt)
        
        coef.string = gsub(" ",".",gsub("-","over",coef))
        tt$Comparison = coef.string

        #add gene annotation
        tt.anno = merge(tt, rat.anno, by.x = "Gene", by.y = "external_gene_name", all.x = T)

        res.df = rbind(res.df, tt.anno)
    }

    column.order = c("Comparison","Gene","logFC","AveExpr","t","P.Value","adj.P.Val")
    res.df = res.df[,column.order]

    #rename columns
    names(res.df) = c("Comparison","Gene","Log2FoldChange","Avg.Expr","T.statistic","Nominal.PValue","BH.PValue")

    return(res.df)
}




#DESeq2 wrapper
de.deseq2 = function(
    counts = NULL, #pseudobulk count matrix
    samples, #sample phenotype metadata
    res.dir, #output location
    cluster = NULL #cluster label for output
){

    #create output folder
    res.file.base = paste0(res.dir,"/DESeq2/",cluster)
    if(!dir.exists(res.file.base)){
        dir.create(res.file.base, recursive = T)
    }
  
    #create DESeq2 object
    dds = DESeqDataSetFromMatrix(countData = round(counts), colData = samples, design = ~ 0+Status)
  
    #perform de testing
    dds = DESeq(dds)
    
    #generate regularized log counts
    rlg = rlog(dds,blind=T)
    lognorm = assay(rlg)
    colnames(lognorm) = paste0(colnames(lognorm), ".RLog")
    write.csv(lognorm, file = paste0(res.file.base,"/","Cluster",cluster,".RlogCounts.csv"))
  
    #create rlog counts scaled across samples for each gene - used for plotting
    scaled.rlog.cts.file = paste0(res.file.base,"/","Cluster.",cluster,".RlogCounts.Scaled.csv")
    cts.scaled = t(scale(t(lognorm)))
    write.csv(cts.scaled, file = scaled.rlog.cts.file)
  
    #combine results from each contrast
    res.df = data.frame()
    for(coef in coefs){
        B = strsplit(coef," - ")[[1]][1]
        A = strsplit(coef," - ")[[1]][2]
        print(coef)
    
        tt <- data.frame(
            results(dds, 
                    contrast=c("Status",B,A),
                    alpha=0.05,
                    pAdjustMethod="BH"
            )
        )

        #baseMean log2FoldChange     lfcSE       stat    pvalue      padj
        names(tt) = c("AveExpr","Log2FoldChange","Log2FC.StdError","Statistic","Nominal.PValue","BH.PValue")  
        tt$Gene = rownames(tt)
        coef.string = gsub(" ",".",gsub("-","over",coef))
        tt$Comparison = coef.string
    
        #add gene annotation
        tt.anno = merge(tt, rat.anno, by.x = "Gene", by.y = "external_gene_name", all.x = T)
    
        res.df = rbind(res.df, tt.anno)
    }
    return(res.df)
}





#load labeled data
so.list = readRDS("../RData/2HC,NC,CDH.labeled.list.rds")


#perivascular fibroblasts has 2 clusters in 2HC which should be combined
xtabs(~Cell.Type.1+sample,data=so.list[["2HC"]]@meta.data)

so.list[["2HC"]]@meta.data[so.list[["2HC"]]$Cell.Type.1 == "Perivascular Fibroblasts Additional","Cell.Type.1"] = "Perivascular Fibroblasts"

xtabs(~Cell.Type.1+sample,data=so.list[["2HC"]]@meta.data)



#calculate at both granular and global cell type assignments
idents = c("Cell.Type.1","Cell.Type.2")

#list to contain pseudobulk counts for each sample and cluster
pseudobulk.counts = list()


for(ident in idents){
    print(ident)
    Idents(so) = ident

    #iterate over clusters 
    for(cl in unique(so@meta.data[,ident])){
        print(cl)
    
        #subset cluster
        so.cl = subset(so, idents = cl)

        #grab raw counts
        cl.cts = GetAssayData(so.cl, assay = "RNA", slot = "counts")

        #grab phenotype metadata
        pheno.cols = c("sample","Donor")

        pheno = so.cl@meta.data[,pheno.cols]

        #skip cluster unless found in all 4 samples
        if(length(unique(md$sample)) == 4){

            #check counts
            xt = xtabs(~sample+Donor,data=pheno)
            print(xt)
            
            #aggregate
            pseudobulk = aggregate.Matrix(t(cl.cts), groupings = pheno, fun = "sum")    
            print(dim(pseudobulk))      

            #transpose to get genes as rows
            pseudobulk = t(as.matrix(pseudobulk))
            dim(pseudobulk)
            
            pseudobulk.counts[[ident]][[cl]] = pseudobulk
        }
    }
}


#apply DGE
for(ident in idents){
    print(ident)
  
    #output folder
    res.dir = paste0(comparisons.dir,"/",ident)

    #loop across clusters
    cl.res = foreach(cl = names(pseudobulk.counts[[ident]])) %dopar% {
        print(cl)
    
        #count input
        counts = pseudobulk.counts[[ident]][[cl]]

        #phenotype input
        samples = data.frame("Sample"=colnames(counts))
    
        #parse metadata
        samples$Sample = gsub("_E215","",samples$Sample)
        samples$Status = sapply(strsplit(samples$Sample,"_"),"[",2L)
        samples$Donor = sapply(strsplit(samples$Sample,"_"),"[",3L)  

        #update experimental group to simpler names 
        samples$Status = ifelse(samples$Status %in% c("HealthyControl"), "HC", samples$Status)
        samples$Status = ifelse(samples$Status %in% c("NitrofenControl"), "NC", samples$Status)
        samples$Status = ifelse(samples$Status %in% c("NitrofenCDH"), "CDH", samples$Status)
        
        rownames(samples) = samples$Sample
        samples = samples[order(samples$Status),]
    
        #create colors for MDS plot
        statuses = as.vector(unique(samples$Status))
        colors = list()
        colors[[statuses[1]]] = "red"
        colors[[statuses[2]]] = "green"
        colors[[statuses[3]]] = "blue"
        
        #expand colors to complete samples metadata table
        for(status in names(colors)){

            #get number of times color should be repeated    
            status.nsamples = nrow(subset(samples, Status == status))
            
            #add color repeats
            colors[[status]] = rep(colors[[status]], status.nsamples)
        }
        samples$Color = as.vector(unlist(colors))
    
        #grab counts for samples in metadata
        counts = counts[,rownames(samples)]
    
        #run limma DGE
        limma.de.res = de.limma(
            counts = counts,
            samples = samples,
            res.dir = res.dir,
            cluster = gsub("\\/","_",cl) #some cell type labels have forward slashed, so replace with underscores
        )
        
        deseq2.de.res = de.deseq2(
            counts = counts,
            samples = samples,
            res.dir = res.dir,
            cluster = gsub("\\/","_",cl) #some cell type labels have forward slashed, so replace with underscores
        )

        limma.de.res$Cluster = cl
        deseq2.de.res$Cluster = cl

        return(
            list(
                "Limma"=limma.de.res,
                "DESeq2"=deseq2.de.res
            )
        )
    }
    names(cl.res) = names(pseudobulk.counts[[ident]])


    #save genelists
      

    #need to invert results before concatenating
    results = list()
    for(cl in names(cl.res)) {
        print(cl)

        #iterate over methods
        for(method in c("Limma","DESeq2")){
    
            #invert results
            results[[method]][[cl]] = cl.res[[cl]][[method]]
        }
    }

    for(method in names(results)){
        
        #concatenate clusters' results
        method.res = list.rbind(results[[method]])

        #re-order columns in output

        #vector of column positions
        column.indices = c(1:ncol(method.res))

        #cluster column position
        cl.index = grep("Cluster",colnames(method.res))

  
        #make cluster column first
        first.column = c(cl.index)
        other.column.indices = column.indices[!column.indices %in% first.column]
        new.column.indices = c(first.column, other.column.indices)
        method.res =  method.res[,new.column.indices]

        #save output
        output.filename = paste0(method,".2HC,NC,CDH.labeled.Pseudobulk.Donors.",ident,".GeneLists.csv")
    
        #remove rows with NA values for adjusted p-value
        method.res.nona = subset(method.res, !is.na(BH.PValue))
    
        write.csv(method.res.nona, file = paste0(res.dir,"/",method,"/",output.filename), row.names = F)
    
    }
}

























