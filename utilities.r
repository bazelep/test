
#####################################################################store.counts Function
store.counts = function(
    so = NULL, #seurat object
    stage = "Initial", #qc stage
    ID.col = "sample"
) {

    counts = c()

    DefaultAssay(so) = "RNA"
    Idents(so) = ID.col
    samples = names(table(Idents(so)))
 
    for(s in samples){ 
        print(s)
  
        s.so = subset(so, idents = s)
  
        n.genes = dim(s.so)[1]
        n.cells = dim(s.so)[2]

        #Genes
        mean.genes.per.cell = mean(s.so@meta.data$nFeature_RNA)
        mean.percentage.genes.per.cell = mean(s.so@meta.data$nFeature_RNA/n.genes) 
    
        #UMI
        mean.umi.per.cell = mean(s.so@meta.data$nCount_RNA) 
    
        #MT
        mean.mito.ratio = mean(s.so@meta.data$mitoRatio)

        #Ribo
        mean.ribo.ratio = mean(s.so@meta.data$riboRatio)
    
        #Complexity
        mean.complexity = mean(s.so$log10GenesPerUMI)
  
        donor.cts = "NA"
        doublet.cts = "NA"

        if("Donor" %in% names(s.so@meta.data)){
            donor.tbl = table(s.so@meta.data$Donor)
            print(donor.tbl)
              donor.cts = paste0(donor.tbl, collapse=":")
        }

        if("Doublet" %in% names(s.so@meta.data)){
            doublet.cts = nrow(subset(s.so@meta.data, Doublet == 1))
        }

        counts = rbind(
            counts,
            c(
            "Stage" = stage,
            "Sample" = s,
            "Number.Cells" = n.cells,
            "Number.Genes" = n.genes,
            "Mean.Detected.GenesPerCell" = mean.genes.per.cell,
            "Mean.Proportion.Detected.GenesPerCell" = mean.percentage.genes.per.cell,
            "Mean.UMIPerCell" = mean.umi.per.cell,
            "Mean.Mitochondrial" = mean.mito.ratio,
            "Mean.Ribosomal" = mean.ribo.ratio,
            "Mean.Complexity" = mean.complexity,
            "Donor.Counts" = donor.cts,
            "Doublet.Counts" = doublet.cts
            )
        )
    }
    
    return(data.frame(counts))
}
#####################################################################end store.counts function






#############################################################qc.plots function
qc.plots = function(so, #seurat object
                    stage = "Initial", #QC stage 
                    
                    f.genes = list( #gene x-intercept and limits 
                      x.min = NULL,
                      x.max = NULL,
                      x.intercept = 250
                    ), 
                    f.umi = list( #umi x-intercept and limits
                      x.min = NULL,
                      x.max = NULL,
                      x.intercept = 500
                    ), 
                    f.mt = 0.2, #MT proportion x-intercept
                    f.rb = 0.2, #Ribo proportion x-intercept
                    f.complexity = 0.8, #Complexity x-intercept

                    legend = T,
                    qc.dir = "../Graphs/QC",
                    ID.col = "sample"
) {
  

  if(!dir.exists(qc.dir)){
      dir.create(qc.dir)
  }

  if(legend){
    legend.string = ""
  } else {
    legend.string = "+NoLegend()"
  }
  
  metadata = so@meta.data
  
  DefaultAssay(so) = "RNA"

  if(is.null(f.genes$x.min)){
    f.genes$x.min = 1
  }
  if(is.null(f.umi$x.min)){
    f.umi$x.min = 1
  }

  if(is.null(f.genes$x.max)){
    f.genes$x.max = as.numeric(max(metadata$nFeature_RNA))
  }
  if(is.null(f.umi$x.max)){
    f.umi$x.max = as.numeric(max(metadata$nCount_RNA))
  }


  # Visualize the number UMIs/transcripts per cell
  pdf(file=paste0(qc.dir,"/UMIs Per Cell - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(color=eval(parse(text=ID.col)), x=nCount_RNA, fill= eval(parse(text=ID.col)))) + 
    coord_cartesian(xlim = c(f.umi$x.min, f.umi$x.max)) +
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = f.umi$x.intercept)
  if(!legend){p = p + NoLegend()}
  print(p)
  dev.off()


  # Visualize the distribution of genes detected per cell via histogram
  pdf(file=paste0(qc.dir,"/Genes Per Cell - Histogram - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(color=eval(parse(text=ID.col)), x=nFeature_RNA, fill= eval(parse(text=ID.col)))) + 
    coord_cartesian(xlim = c(f.genes$x.min, f.genes$x.max)) +
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = f.genes$x.intercept)
  if(!legend){p = p + NoLegend()}
  print(p)
  dev.off()


  # Visualize the distribution of genes detected per cell via boxplot
  pdf(file=paste0(qc.dir,"/Genes Per Cell - Boxplot - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(x=eval(parse(text=ID.col)), y=log10(nFeature_RNA), fill=eval(parse(text=ID.col)))) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells vs NGenes")
  if(!legend){p = p + NoLegend()}
  print(p)
  dev.off()


  # Visualize the number of cell counts per sample
  pdf(file=paste0(qc.dir,"/Cells Per Sample - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(x=eval(parse(text=ID.col)), fill=sample)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  if(!legend){p = p + NoLegend()}
  print(p)
  dev.off()

  # Visualize the correlation between genes detected and number of UMIs and determine 
  # whether strong presence of cells with low numbers of genes/UMIs
  pdf(file=paste0(qc.dir,"/Genes Per UMI - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
    geom_point(size=0.3) + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() + 
    geom_vline(xintercept = f.umi$x.intercept) +
    geom_hline(yintercept = f.genes$x.intercept) +
    facet_wrap(~eval(parse(text=ID.col)))
  print(p)
  dev.off()

  # Visualize the distribution of mitochondrial gene expression detected per cell
  pdf(file=paste0(qc.dir,"/MT Ratio - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(color=eval(parse(text=ID.col)), x=mitoRatio, fill=eval(parse(text=ID.col)))) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = f.mt)
  if(!legend){p = p + NoLegend()}
  print(p)
  dev.off()

  # Visualize the distribution of ribosomal gene expression detected per cell
  pdf(file=paste0(qc.dir,"/Ribosomal Ratio - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(color=eval(parse(text=ID.col)), x=riboRatio, fill=eval(parse(text=ID.col)))) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = f.rb)
  if(!legend){p = p + NoLegend()}
  print(p)
  dev.off()

  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  pdf(file=paste0(qc.dir,"/Complexity - ",stage,".pdf"))
  p = ggplot(metadata,
    aes(x=log10GenesPerUMI, color = eval(parse(text=ID.col)), fill=eval(parse(text=ID.col)))) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = f.complexity)
  if(!legend){p = p + NoLegend()}
  print(p)
  dev.off()  
}
#############################################################end qc.plots function


#####################################################################SCT.normalize function
SCT.normalize = function(
    so,
    cell.cycle.genes = NULL, 
    MT = TRUE,
    variable.features.n = 3000,
    n.cells = 5000
) {

    so.list = SplitObject(so, split.by = "sample")

    if(MT){
        vars = c("mitoRatio")
    } else {
        vars = NULL
    }

    for (i in 1:length(so.list)) {
        print(names(so.list)[i])

        #needed for cell cycle scoring
        so.list[[i]] = NormalizeData(so.list[[i]], verbose = TRUE)

        if(!is.null(cell.cycle.genes)){
            print("...Found cell cycle genes variable...")
            vars = c("mitoRatio","S.Score", "G2M.Score")

            so.list[[i]] = CellCycleScoring(so.list[[i]], 
                                    g2m.features=cell.cycle.genes$g2m_genes, 
                                    s.features=cell.cycle.genes$s_genes)
        }

        so.list[[i]] = SCTransform(so.list[[i]], 
                            vars.to.regress = vars, 
                            variable.features.n = variable.features.n,
                            ncells = n.cells,
                            method = "qpoisson"
                            )

        #needed for reciprocal PCA
        so.list[[i]] = RunPCA(so.list[[i]], assay = "SCT")
    }
    
    return(so.list)
}
#####################################################################end SCT.normalize function




#####################################################################get.annotations function
get.annotations = function(
  species = "Homo_sapiens", 
  snapshot.date = NULL
){

  species = gsub("_"," ",species)
  
  #connect to AnnotationHub
  ah = AnnotationHub()

  snapshotDate(ah) = snapshot

  #access the Ensembl database for organism
  ahDb = query(
    ah,
    pattern = c(species, "EnsDb"), 
    ignore.case = TRUE
  )
  id = rownames(tail(mcols(ahDb),n=1))

  #download the appropriate Ensembldb database
  edb = ah[[id]]
  
  #extract gene-level information from database
  annotations = genes(edb, return.type = "data.frame")
  
  #select annotations of interest
  annotations = annotations[,c("gene_id","gene_name","seq_name","gene_biotype","description")]

  return(annotations)
}
#####################################################################end get.annotations function


#####################################################################make.clusters function
make.clusters = function(
  so, #seurat object
  resolutions, #vector of resolutions for Leiden clustering
  sample.set = "ALL", #tag for output labeling
  split.by = NULL, #DimPlot() split.by; use this or grid.samples option
  grid.samples = NULL, #which samples to show in plot grid
  pca.dims = 30,
  assay = "integrated",
  reductions = NULL, #tsne or umap
  reduction.dims = c(1,2), #how many tsne/umap dimensions to plot
  output.dir.name = "Clusters",
  sample.col = "sample" #sample metadata column to use
) {

  cluster.plots.dir = paste0("../Graphs/",output.dir.name)
  cluster.counts.dir = paste0("../OutputData/",output.dir.name)

  if(!dir.exists(cluster.plots.dir)){
     dir.create(cluster.plots.dir)
  }

  if(!dir.exists(cluster.counts.dir)){
     dir.create(cluster.counts.dir)
  }

  if(!"pca" %in% names(so@reductions)){
    DefaultAssay(so) = assay
    so = RunPCA(so)
    so = RunUMAP(so, dims = 1:pca.dims,reduction = "pca", return.model = T)
    so = RunTSNE(so)
  }

  if(is.null(reductions)){
    reductions = c("tsne","umap")
  }

  #determine the K-nearest neighbor graph
  DefaultAssay(so) = assay
  so = FindNeighbors(so, dims = 1:pca.dims)

  #find clusters in graph
  so = FindClusters(so, resolution = resolutions)


  #store cluster cell counts and plot tsne/umap for each clustering resolution
  for(res in resolutions){
  
    ID.col = paste0(assay,"_snn_res.",res)

    Idents(so) = ID.col
  
    #save cluster cell counts
    cell.counts = FetchData(so,vars = c("ident", sample.col))

    cell.counts = as.matrix(xtabs(~ eval(parse(text=sample.col)) + ident,data=cell.counts))
        
    cell.counts.filename = paste0(
      cluster.counts.dir,"/",
      sample.set,
      " - Resolution ",res,
      " - Cluster Cell Counts by Sample.csv"
    )
    
    write.csv(cell.counts, file = cell.counts.filename)


    #generate plots, either as a grid of samples or using split.by option in DimPlot(), or don't split if NULL
    plot.filename.base = paste0(
      cluster.plots.dir,"/", 
      sample.set, 
      " - Resolution ",res,
      " - Dims",paste(reduction.dims,collapse=",")
    )
    
    if(!is.null(grid.samples)){

      Idents(so) = sample.col

      for(red in reductions){
        grid.plots = list()

        for(sam in grid.samples){
          print(sam)
          sam.so = subset(so, idents = sam)

          Idents(sam.so) = paste0(assay,"_snn_res.",res)

          grid.plots[[sam]] = DimPlot(sam.so, 
                            label = T,
                            label.size = 2, 
                            reduction = red,
                            dims = reduction.dims
                            ) + NoLegend() + ggtitle(sam)
            
        }
        pdf(file=paste0(plot.filename.base, toupper(red),".pdf"))
        p = DimPlot(so,
              reduction = red,
              label = TRUE,
              group.by = ID.col,
              label.size = 6,
              dims = reduction.dims) + NoLegend()
        print(p)
        p = DimPlot(so,
              reduction = red,
              label = F,
              group.by = sample.col,
              label.size = 6,
              dims = reduction.dims)
        print(p)
        p = wrap_plots(grid.plots)
        plot(p)
        dev.off()
      }
    } else {
      for(red in c("tsne","umap")){
        pdf(file=paste0(plot.filename.base, toupper(red),".pdf"))
        p = DimPlot(so,
              reduction = red,
              label = TRUE,
              split.by = split.by,
              label.size = 6,
              dims = reduction.dims) + NoLegend()
        print(p)
        dev.off()
      }
    }
  }
  return(so)
}
#####################################################################end make.clusters function




#####################################################################get.conserved function
get.conserved = function(
  cluster, #cluster id
  so, #seurat object
  annotations, #annotations dataframe to add to marker lists
  n.groups, #number of groups expected in data - to test for missing
  direction = "positive", #fold-change direction to filter by
  group.samples.by = "sample"
){
    
    results = tryCatch({
      
      FindConservedMarkers(so,
                           ident.1 = cluster,
                           only.pos = ifelse(direction == "both", F, T),
                           min.cells.feature = 1,
                           min.cells.group = 1,
                           logfc.threshold = 0.25,
                           test.use = "wilcox",
                           min.pct = 0.1,
                           grouping.var = group.samples.by
                           ) 
    }, error = function(err){
      print(paste("FCM ERROR:",err))  
    })
    
    print(paste0("...finished FCM for cluster ",cluster,"..."))

    #check if FCM failed 
    if( (class(results) == 'data.frame') && (nrow(results) > 0) ){
      print("..valid results...")
      
      #add annotations to results
      results$gene = rownames(results)
      results = merge(
        results, 
        unique(annotations[,c("gene_name","description")]),
        by.x = "gene", by.y = "gene_name"
      )
      results$cluster_id = cluster

      #check that expected number of groups is equal to returned
      # number of groups, in case of low cell count exclusion
      results.colnames = names(results)
      groups.included = results.colnames[grep("avg_log2FC",results.colnames)]
      
      if(length(groups.included) == n.groups){
        print("...group lengths match...")

        if (n.groups > 1){
          print("...more than one group...")
          print(paste0("Number of groups (",length(groups.included), ") with conserved cluster ",cluster,
                      " equal to number of expected groups (",n.groups,")..."))

          print(head(results))
          print(groups.included)

          #calculate average fold-change
          fc = as.data.frame(results[,groups.included])
          names(fc) = groups.included

          fc.avg = function(row){
            #row = 1
            row.sum = sum(fc[row,])
            row.sum / ncol(fc)
          }
          results$avgFC = sapply(1:nrow(fc), fc.avg)
          
        } else {
          print("..only 1 sample...")
          fc.col = grep("avg_log2FC",results.colnames,value=T)
          pv.col = grep("p_val_adj",results.colnames,value=T)

          results$avgFC = results[,fc.col]
          results$minimump_p_val = res[,pv.col]
        }

        return(results)
      } else { #number of expected groups differs from returned groups
        print(
          paste0(
            "Number of groups (",length(groups.included), ") with conserved cluster ",cluster,
            " not equal to number of samples (",n.groups,")..."
          )
        )

        return(NULL)
      }
    } else { #results not dataframe or has 0 rows
      print("..invalid data:")
      print(results)
      return(NULL)
    }
  }
#####################################################################end get.conserved function


#####################################################################get.markers function
get.markers <- function(
  so, #seurat object
  annotations, #annotations to add to marker list
  sample.set = "ALL", #tag for output
  direction = "positive", #fold-change direction to include
  ID.col = NULL,
  group.samples.by = "sample", #cell grouping metadata column
  output.dir = "../OutputData/MarkerGenes/"
){

  if(!dir.exists(output.dir)){
    dir.create(output.dir)
  }

  #track how many sample groups present in data, in case some are excluded
  # for having too few cells in a given cluster
  Idents(so) = group.samples.by
  n.groups = length(unique(Idents(so)))

  if(!is.null(ID.col)){
    Idents(so) = ID.col
  } else {
    stop("...ID.col missing")
  }

  print("...cluster cell counts:")
  print(table(Idents(so)))

  #identities = as.character(sort(as.integer(unique(as.character(FetchData(so, vars = c("ident"))$ident)))))
  #clusters = as.character(sort(as.integer(unique(Idents(so)))))
  clusters = unique(Idents(so))
  print("...unique clusters:")
  print(as.character(clusters))
  
  if(length(clusters) > 1){

    so$cluster_id = Idents(so)

    print("...sample group by cluster counts:")
    xt = xtabs(
        ~ eval(parse(text=group.samples.by)) + cluster_id, 
        data = so@meta.data)
    
    print(xt)

    #filter out clusters not found in all sample groups (not conserved)
    clusters.keep = clusters[colSums(xt > 1) == n.groups]
    
    print("...clusters found in all sample groups:")
    print(clusters.keep)

    if(length(clusters.keep) > 1){

      #gather conserved markers for each cluster
      conserved.markers = data.frame()

      for(cl in clusters.keep){
        print(paste0("...looking for conserved marker genes in cluster: ",cl))
       
        cl.cm = get.conserved(
          cluster = cl, 
          so = so, 
          annotations = annotations, 
          n.groups = n.groups, 
          direction = direction,
          group.samples.by = group.samples.by
        )

        if( (!is.null(cl.cm)) & (class(cl.cm) == "data.frame") ){
          print(head(cl.cm))
          conserved.markers = rbind(conserved.markers, cl.cm)
        }
      }


      #save marker list
      markers.file = paste0(output.dir,"/",sample.set, "-Markers-",ID.col,".csv")
      write.csv(conserved.markers, file = markers.file,row.names=F)
      
      #create an aggregated version of markers, where each column
      # represents a cluster and contains the minimum p-value and averaged fold-change
      agg = conserved.markers[,c("cluster_id","gene","description","avgFC","minimump_p_val")]
      
      #cleanup description
      agg$description = sub(" \\[Source.*$", "", agg$description, perl = T)
      
      #create combined result column
      agg$result = paste0(agg$gene," (",agg$description,")(",round(agg$avgFC,digits=2),")(",agg$minimump_p_val,")")
      
      #subset columns
      agg = agg[,c("cluster_id","result","avgFC","minimump_p_val")]

      #since each cluster may have differing numbers of genes (rows), need
      # to populate "extra" rows in order to cbind correctly
      agg.cluster.ids = sort(unique(agg$cluster_id))

      #get number of genes (rows) for each cluster
      get.cluster.nrow = function(cl){
        nrow(subset(agg, cluster_id == cl))
      }
      n.hits = sapply(agg.cluster.ids, get.cluster.nrow)

      #get maximum number of genes (rows)
      max.n.hits = max(n.hits)
         
      #sort results by average fold-change and then minimum p-value
      # and extract 'result' column
      sort.results = function(cl){
        cl.hits = subset(agg, cluster_id == cl)
    
        cl.hits.results = cl.hits[
          order(
            cl.hits$avgFC,
            cl.hits$minimump_p_val, 
            decreasing = c(T,F)
          )
        ,]$result

        #add NAs to end if shorter than max.n.hits
        length(cl.hits.results) = max.n.hits
        return(cl.hits.results)
      }
      
      #sort results and extract aggregated column
      agg.results = list.cbind(lapply(agg.cluster.ids, sort.results))

      #label each column with cluster id
      colnames(agg.results) = paste0("Cluster_",agg.cluster.ids)

      write.csv(agg.results, 
          file = paste0(output.dir,"/",sample.set,"-Markers-",ID.col,"-Grouped.csv"),
          row.names=F)

    } else{
      print("...no clusters found in all samples, skipping...")
    }
  } else { #clusters > 1
    print("...<= 1 clusters found, skipping...")
  }
}
#####################################################################end get.markers function

#####################################################################run.liana function
run.liana = function(
  sce = NULL, 
  comparison = NULL, 
  sample = NULL,
  output.dir = NULL
){
    

  liana.res = liana_wrap(sce = sce,
                        idents_col = "Ligand.Receptor.Group"
                        )
  liana.agg = liana_aggregate(liana.res)

  liana.results = list("separate"=liana.res,"aggregate"=liana.agg)
    
  saveRDS(liana.results, file = paste0(output.dir,"/",comparison,".",sample,".results.rds"))
  
}
#####################################################################end run.liana function
