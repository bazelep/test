#Summary: Call germline variants for per-cell SAM files
#Conda: cdh


library(doFuture)

registerDoFuture()

set.seed(100)
options(width=100)




donors.dir = "../OutputData/Donors/"
genome.fasta.file = "../OutputData/Cellranger/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"


#genome fasta needs to be indexed for freebayes
if(!file.exists(paste0(genome.fasta.file,".fai"))){
    cat(paste0("...indexing genome fasta..."))
    system(
        paste0(
            "samtools faidx ",genome.fasta.file
        )
    )
}

barcode.files = list.files(path = donors.dir, pattern = "\\.barcodes\\.csv", full.names = T)


#plan("multisession", workers = 100)






#parallize running scripts

#number of tasks to run in parallel
batch.size = 150 

#store each type of task as a list entry
# - tasks are a vector of command strings for the given list entry
cmds.list = list()


files.res = foreach(barcode.file = barcode.files) %do% {
    print(barcode.file)

    sample.barcodes = read.table(barcode.file, header = F,sep = "\t")$V1

    sample.name = gsub("\\.barcodes\\.csv", "", basename(barcode.file))
    sample.dir = paste0(donors.dir,"/",sample.name)

    sample.barcodes.dir = paste0(sample.dir,"/barcodes")

    sample.sam = paste0(sample.name,".sam")
    sample.sam.header = paste0(sample.sam,".header")

    #overall task (barcode) counter for sample
    counter = 0

    #keep looping until read last sample barcode
    while(counter < length(sample.barcodes)){

        counter = counter + 1

        barcode = sample.barcodes[counter]
    
        #add wait command at end of each batch; if within-batch, add ampersand
        if( 
            ((counter %% batch.size) == 0) 
            | 
            (counter == length(sample.barcodes))
        ){
           tail = "; wait;"
        } else {
            tail = " &";
        }

        barcode.sam = paste0(sample.barcodes.dir, "/", sample.name,".",barcode,".sam")
        barcode.bam = paste0(barcode.sam,".bam")
        barcode.vcf = paste0(barcode.bam,".fb.vcf")
        barcode.vcf.gz = paste0(barcode.vcf,".gz")

        #save each task type to a separate vector
        cmds.list[["make_bam"]] = c(cmds.list[["make_bam"]] , paste0("samtools view -hb ",barcode.sam," > ",barcode.bam,tail))

        cmds.list[["freebayes"]] = c(cmds.list[["freebayes"]] , paste0("freebayes -f ",genome.fasta.file," ",barcode.bam," > ",barcode.vcf,tail))
        
        cmds.list[["update_vcf"]] = c(cmds.list[["update_vcf"]] , paste0("perl ./update_vcf_sample_name.pl ",barcode.vcf,tail))
        
        cmds.list[["mv_update"]] = c(cmds.list[["mv_update"]] , paste0("mv ",barcode.vcf,".updated ",barcode.vcf,tail))
        
        cmds.list[["bgzip"]] = c(cmds.list[["bgzip"]] , paste0("bgzip -f ",barcode.vcf,tail))
        
        cmds.list[["index"]] = c(cmds.list[["index"]] , paste0("bcftools index -f ",barcode.vcf.gz,tail))
    
    }
}



cmds.res = foreach(cmd.type = names(cmds.list)) %do% {
    print(cmd.type)

    #very last command should have 'wait' after it
    cmds.list[[cmd.type]][[length(cmds.list[[cmd.type]])]] = gsub("&",";wait",cmds.list[[cmd.type]][[length(cmds.list[[cmd.type]])]])

    cmds.file = paste0("../OutputData/Donors/",cmd.type,".sh")
    write.table(cmds.list[[cmd.type]], file = cmds.file, row.names=F,col.names=F,quote=F)

    system(paste0("sh ",cmds.file))
}




