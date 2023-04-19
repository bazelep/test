#Summmary: Prepare Cellranger binary and reference files
#Conda: cdh



#create folder for results
cellranger.dir = "../OutputData/Cellranger/"
cellranger.binary = paste0(cellranger.dir, "/cellranger-4.0.0/bin/cellranger")

if(!dir.exists(cellranger.dir)) {dir.create(cellranger.dir,recursive=T)}


#download Cellranger binary
system(
    paste0(
        'wget ',
        '-o ',cellranger.dir,'/cellranger.wget.log',
        '-O ',cellranger.dir,'/cellranger-4.0.0.tar.gz ',
        '"https://cf.10xgenomics.com/releases/cell-exp/cellranger-4.0.0.tar.gz?Expires=1680849401&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci00LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODA4NDk0MDF9fX1dfQ__&Signature=UIkKoxACcldjSttXGC-vhq8dVcJnr7zNFyhicv95x8DhaXLI9Azmaskre05AM7wVi6BHQ7kyf-V-0l1v9FLREgkgxqO9clEfVoFG690DW2GwdUdF5Hhi5Ao6aS1lP9j5QiLGSUaXe19YsHG1KNd~f4sB1~FpGhpfwl4hf2SDG3YoXAIS0ajftBr-9RoJ1iQhHwL8o7MZcoxHSxjYldXhrEuWsUO0TJXPXaninshZ0mjFAPH1-6XZeSNKgeic-oV7-5sybF5D2oe~4bk4nxwimZe-E3zbpCg5xqHIcbEVGmGUYVrzhXffVyHExlpPG-qyWSZX3xrX60Imjqk04e5n9Q__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"'
    )
)
system(
    paste0(
        "tar xzf ",cellranger.dir,"/cellranger-4.0.0.tar.gz -C ",cellranger.dir
    )
)



##Need to create custom reference since only human and mouse available on 10x website


#download genome fasta and gtf
genome.fasta.file.link = "http://ftp.ensembl.org/pub/release-101/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
genome.gtf.file.link = "http://ftp.ensembl.org/pub/release-101/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.101.gtf.gz"
system(
    paste0(
        "wget ",
        "-o ",cellranger.dir,"/fasta.wget.log",
        " -P ",cellranger.dir," "
        genome.fasta.file.link
    )
)
system(
    paste0(
        "wget ",
        "-o ",cellranger.dir,"/gtf.wget.log",
        " -P ",cellranger.dir," "
        genome.gtf.file.link
    )
)
system(
    paste0(
        "gunzip -d ",cellranger.dir,"/*.gz"
    )
)



#Create custom reference
genome.base = paste0(cellranger.dir, "/Rattus_norvegicus.Rnor_6.0")
genome.gtf = paste0(genome.base,".101.gtf")
genome.fasta = paste0(genome.base,".dna.toplevel.fa")

mkgtf.cmd = paste0(
    cellranger.binary," ",
    "mkgtf ",genome.gtf," ", gsub("gtf","",genome.gtf),"filtered.gtf ",
    "--attribute=gene_biotype:protein_coding ",
    "--attribute=gene_biotype:lincRNA ",
    "--attribute=gene_biotype:antisense ",
    "--attribute=gene_biotype:IG_LV_gene ",
    "--attribute=gene_biotype:IG_V_gene ",
    "--attribute=gene_biotype:IG_V_pseudogene ",
    "--attribute=gene_biotype:IG_D_gene ",
    "--attribute=gene_biotype:IG_J_gene ",
    "--attribute=gene_biotype:IG_J_pseudogene ",
    "--attribute=gene_biotype:IG_C_gene ",
    "--attribute=gene_biotype:IG_C_pseudogene ",
    "--attribute=gene_biotype:TR_V_gene ",
    "--attribute=gene_biotype:TR_V_pseudogene ",
    "--attribute=gene_biotype:TR_D_gene ",
    "--attribute=gene_biotype:TR_J_gene ",
    "--attribute=gene_biotype:TR_J_pseudogene ",
    "--attribute=gene_biotype:TR_C_gene",
    " &> ", genome.gtf,".mkgtf.log"
)
mkgtf.cmd

system(mkgtf.cmd)

genome.output.dir = paste0("refdata-gex-",basename(genome.base),".101")

mkref.cmd = paste0(
    cellranger.binary," ",
    "mkref ",
    "--genome=",genome.output.file, " ",
    "--fasta=",genome.fasta," ",
    "--genes=",gsub("gtf","",genome.gtf),"filtered.gtf ",
    "--memgb=200 ",
    "--nthreads=40 ",
    "&> ", genome.base,".mkref.log"
)
mkref.cmd

system(mkref.cmd)

#coudn't specify output path so move
system(
    paste0(
        "mv ",genome.output.dir," ",cellranger.dir,";",
        "mv Log.out ",cellranger.dir
    )
)






