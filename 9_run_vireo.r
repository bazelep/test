#Summary: Run Vireo on per-cell germline variants
#Conda: cdh


library(doFuture)

registerDoFuture()

set.seed(100)
options(width=100)




donors.dir = "../OutputData/Donors/"
n.donors = 5 #all samples have 5 donors



barcode.files = list.files(path = donors.dir, pattern = "csv", full.names = T)



files.res = foreach(barcode.file = barcode.files) %do% {
    print(barcode.file)

    sample.name = gsub("\\.barcodes\\.csv", "", basename(barcode.file))
    sample.dir = paste0(donors.dir,"/",sample.name)
    sample.vireo.dir = paste0(sample.dir,"/vireo.out")
    
    if(!dir.exists(sample.vireo.dir)){
        dir.create(sample.vireo.dir)
    }

    sample.barcodes.dir = paste0(sample.dir,"/barcodes")

    #merge all the vcfs
    cmd = paste0(
        "ls ",sample.barcodes.dir,"/*.vcf.gz > ",sample.dir,"/vcf.list.txt"
    )
    print(cmd)
    system(cmd)

    cmd = paste0(
        "bcftools merge -l ",sample.dir,"/vcf.list.txt -O z -o ",sample.dir,"/",sample.name,".merged.vcf.gz"
    )   
    print(cmd)
    system(cmd)

    #filter merged vcf for variants with >100 alleles (>100 reads across all cells)
    ###INFO=<ID=AN,Number=1,Type=Integer,Description="Number of Alleles in Samples with Coverage">
    cmd = paste0(
        "bcftools view -i 'INFO/AN > 100' ",
        sample.dir,"/",sample.name,".merged.vcf.gz | bgzip -c > ",
        sample.name,".merged.AN100.vcf.gz"
    )
    print(cmd)
    system(cmd)

    #run vireo on filtered vcf
    cmd = paste0(
        "vireo -c ",sample.name,".merged.AN100.vcf.gz -N ",n.donors," -o ",sample.vireo.dir,
        " &> ",sample.vireo.dir,"/",sample.name,".vireo.log"
    )
    print(cmd)
    system(cmd)
}






