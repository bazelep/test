#Summmary: Run Cellranger to align/quantify to generate counts
#Conda: cdh

library(doFuture)

registerDoFuture()
plan("sequential")


cellranger.dir = "../OutputData/Cellranger/"
cellranger.binary = paste0(cellranger.dir,"/cellranger-4.0.0/bin/cellranger")
fastqs.base.dir = "../OutputData/FastQ/"
counts.dir = "../OutputData/Counts/"
genome.base = paste0(cellranger.dir, "/Rattus_norvegicus.Rnor_6.0")
reference.dir = paste0(cellranger.dir,"/refdata-gex-",basename(genome.base),".101")

#Create output counts folder
if(!dir.exists(counts.dir)){dir.create(counts.dir,recursive=T)}

# system(
#     paste0(
#         cellranger.binary," ",
#         "count ",
#         "--help"
#     )
# )

bcl.dirs = list.files(path="../Data/Raw_Data", pattern = "XX", full.names=T)

#Not including E18.5 or E14.5 anymore, but if do, add --force-cells 8000 for E18.5


bcl.res = foreach(bcl.dir = bcl.dirs) %do% {
    print(bcl.dir)

    fastqs.dir = list.files(
        path = paste0(fastqs.base.dir, "/", basename(bcl.dir), "/outs/fastq_path/"), 
        pattern = "XX", 
        full.names=T
    )
    fastqs.dirs = list.files(path=fastqs.dir,full.names=T)

    sample.res = foreach(fastq.dir = fastqs.dirs) %do% {
        print(fastq.dir)

        id = basename(fastq.dir)
        print(id)

        count.cmd = paste0(
            cellranger.binary," ",
            "count ",
            "--id ",id," ",
            "--transcriptome ",reference.dir," ",
            "--localcores 100 ",
            "--localmem 500 ",
            "--fastqs ",fastq.dir," ",
            "&> ",counts.dir,"/",id,".count.log"
        )
        print(count.cmd)

        system(count.cmd)
        system(
            paste0(
                "mv ",id," ",counts.dir
            )
        )
        system(
            paste0(
                "mv __",id,".mro ",counts.dir
            )
        )

    }
}


#newer samples came as fastq

fastq.res = foreach(fastq.dir = c("Sample_FC.JA.2998_1","Sample_FC.JR.3032_2")) %do% {
    print(fastq.dir)

    id = gsub("\\.","_",basename(fastq.dir))
    print(id)

    count.cmd = paste0(
        cellranger.binary," ",
        "count ",
        "--id ",id," ",
        "--transcriptome ",reference.dir," ",
        "--localcores 100 ",
        "--localmem 500 ",
        "--fastqs ",fastqs.base.dir,"/",fastq.dir," ",
        "&> ",counts.dir,"/",id,".count.log"
    )
    print(count.cmd)

    system(count.cmd)
    system(
        paste0(
            "mv ",id," ",counts.dir
        )
    )
    system(
        paste0(
            "mv __",id,".mro ",counts.dir
        )
    )
}


#multi-fastq sample
fastq.dirs =  c("Sample_FC.JR.3032_1_1","Sample_FC.JR.3032_1_2")
fastqs.string = paste(paste0(fastqs.base.dir,"/",fastq.dirs),collapse=",")



id = "Sample_FC_JR_3032_1"
print(id)

count.cmd = paste0(
    cellranger.binary," ",
    "count ",
    "--id ",id," ",
    "--transcriptome ",reference.dir," ",
    "--localcores 100 ",
    "--localmem 500 ",
    "--fastqs ",fastqs.string," ",
    "&> ",counts.dir,"/",id,".count.log"
)
print(count.cmd)

system(count.cmd)

system(
    paste0(
        "mv ",id," ",counts.dir
    )
)
system(
    paste0(
        "mv __",id,".mro ",counts.dir
    )
)



