#Summmary: Run Cellranger to generate fastq files
#Conda: cdh


library(doFuture)

registerDoFuture()
plan("sequential")


cellranger.dir = "../OutputData/Cellranger/"
cellranger.binary = paste0(cellranger.dir,"/cellranger-4.0.0/bin/cellranger")
fastqs.base.dir = "../OutputData/FastQ/"



if(!dir.exists(fastqs.base.dir)){dir.create(fastqs.base.dir,recursive=T)}


# system(
#     paste0(
#         cellranger.binary," ",
#         "mkfastq ",
#         "--help"
#     )
# )



bcl.dirs = list.files(path="../Data/Raw_Data", pattern = "XX", full.names=T)

bcl.res = foreach(bcl.dir = bcl.dirs) %do% {
    print(bcl.dir)

    simplecsv.files = list.files(path=bcl.dir,pattern="simplecsv",full.names=T)

    if(length(simplecsv.files) > 1){
        stop(paste0("...more than one sample sheet:",simplecsv.files))
    } else {
        simplecsv.file = simplecsv.files[1]
    }
    print(simplecsv.file)

    id = basename(bcl.dir)
    print(id)

    mkfastq.cmd = paste0(
        cellranger.binary," ",
        "mkfastq ",
        "--run=",bcl.dir," ",
        "--simple-csv=",simplecsv.file," ",
        "--qc ",
        "--localcores=40 ",
        "--localmem=200 ",
        "--id=",id," ",
        "&> ",fastqs.base.dir,"/",id,".mkfastq.log"
    )
    print(mkfastq.cmd)
    system(mkfastq.cmd)

    system(
        paste0(
            "mv ",id, " ", fastqs.base.dir
        )
    )
    system(
        paste0(
            "mv __",id,".mro ", fastqs.base.dir
        )
    )
}


#newer samples came as fastq

#special cases 

system(
    paste0(
        "cp -r ../Data/Raw_Data/FC.JR.2998/Sample_FC.JA.2998_1/ ../OutputData/FastQ/"
    )
)

system(
    paste0(
        "cp -r ../Data/Raw_Data/FC.JR.3032/Sample_FC.JR.3032_2/ ../OutputData/FastQ/"
    )
)


system(
    paste0(
        "cp -r ../Data/Raw_Data/FC.JR.3032/Sample_FC.JR.3032_1 ../OutputData/FastQ/Sample_FC.JR.3032_1_1"
    )
)

system(
    paste0(
        'mkdir ../OutputData/FastQ/Sample_FC.JR.3032_1_2;',
        'cp -r ../Data/Raw_Data/FC.JR.3032/Sample_FC.JR.3032_1/more\\ reads\\ for\\ Sample_FC.JR.3032_1/* ../OutputData/FastQ/Sample_FC.JR.3032_1_2/'
    )
)



#files came with periods in filenames, which Cellranger doesn't like
fastq.dirs = list.files(path = fastqs.base.dir, pattern = "Sample")

dirs.res = foreach(fastq.dir = fastq.dirs) %do% {
    print(fastq.dir)

    fastq.files = list.files(path = paste0(fastqs.base.dir,"/",fastq.dir))
    print(fastq.files)

    file.res = foreach(fastq.file = fastq.files) %do% {
        print(fastq.file)

        new.fastq.file = paste0(gsub(".fastq.gz","",gsub("\\.","_",fastq.file)),".fastq.gz")

        mv.cmd = paste0(
            "mv ",
            fastqs.base.dir,"/",fastq.dir,"/",fastq.file,
            " ",
            fastqs.base.dir,"/",fastq.dir,"/",new.fastq.file
        )
        print(mv.cmd)
        system(mv.cmd)
    }
}





