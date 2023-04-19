#Summary: Separate per-cell reads from Cellranger BAMs 
#Conda: cdh


library(Seurat)
library(doFuture)

registerDoFuture()

set.seed(100)
options(width=100)

source("utilities.r")



donors.dir = "../OutputData/Donors/"
counts.dir = "../OutputData/Counts/"

#need original sample file names to get BAM files

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


#mapping of original sample file name to barcode file names
sample.map = unique(filtered_seurat@meta.data[,c("orig.ident","sample")])


barcode.files = list.files(path = donors.dir, pattern = "csv", full.names = T)



files.res = foreach(barcode.file = barcode.files) %do% {
    print(barcode.file)

    sample.name = gsub("\\.barcodes\\.csv", "", basename(barcode.file))
    sample.dir = paste0(donors.dir,"/",sample.name)
    sample.barcodes.dir = paste0(sample.dir,"/barcodes")

    if(!dir.exists(sample.barcodes.dir)){
        dir.create(sample.barcodes.dir,recursive = T)
    }

    sample.file.name = subset(sample.map, sample == sample.name)$orig.ident

    sample.bam.dir = paste0(counts.dir, "/", sample.file.name,"/outs")
    sample.bam = paste0(sample.bam.dir,"/possorted_genome_bam.bam")
    sample.sam = paste0(sample.name,".sam")
    sample.sam.header = paste0(sample.sam,".header")


    cmd = paste0(
        "samtools view -H ",
        sample.bam,
        "> ",sample.dir,"/", sample.sam.header
    )
	print(cmd)
    system(cmd)

    cmd = paste0(
        "samtools view -h ",
        sample.bam," ",
        " > ",sample.dir,"/",sample.sam
    )
	print(cmd)
    system(cmd)

    cmd = paste0(
        "time perl ./create_per_cell_sam.pl ",
        sample.name," ",
        barcode.file," ",
        sample.dir,"/",sample.sam," ",
        sample.dir,"/",sample.sam.header," ",
        sample.barcodes.dir," ",
        " &> ",sample.dir,"/create_per_cell_sam.",sample.name,".log"
    )
    print(cmd)
    system(cmd)

}


