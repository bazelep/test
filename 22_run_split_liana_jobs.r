#Summary: Run split up receptor-ligand comparisons with liana.
#Conda: cdh

library(doFuture)

set.seed(100)
options(width=100)


#get full path for run.liana function
project.dir = getwd()

#output folder
ligand.receptor.dir = paste(project.dir,"/../OutputData/LigandReceptor/")
comparisons.dir = paste0(ligand.receptor.dir, "/Comparisons")


#use 'run_liana' function in utilities file
utilities.filename = paste0(getwd(),"/utilities.r")

#get list of SCE objets
sce.files = list.files(path=comparisons.dir,pattern="sce",full.names=T)
sce.files

#iterate of SCE objects (ligand-receptor cell type pairs to be compared) and create individual R scripts
# - to perform liana analysis
sce.res = foreach(sce.file = sce.files) %do% {

    #comparison row number
    comp = strsplit(basename(sce.file),"\\.")[[1]][1]
    
    #status
    status = strsplit(basename(sce.file),"\\.")[[1]][2]

    #code template
    code = paste0(
"
library(SingleCellExperiment)
library(scuttle)
library(liana)
load('",sce.file,"')
sce <- logNormCounts(sce)
source('",utilities.filename,"')
run.liana(sce=sce,comparison='",comp,"',sample='",sample,"',output.dir='",comparisons.dir,"')
"
    ) 

    #save code to file to be run in next script
    write.table(code, file = paste0(comparisons.dir,"/",comp,".",status,".script.R"),row.names=F,col.names=F,quote=F)
}


#parallize running scripts

#maximum comparison row number
max.comps = 465 

#number of jobs to run in parallel
batch.size = 100 

#vector of commands
cmds = c()



for(i in c(0:9)){
    #print(paste0("..i is ",i))

    for(j in c(1:batch.size)){
        #print(paste0("...j is ",j))
    
        #comparison row number 
        comp = (i*batch.size) + j

        #if less than last comparison
        if (comp <= max.comps) {
            #print(paste0("...comp is ",comp))
    
            #iterate over statuses
            for(status in c("2HC","NC","CDH")){
                print(status)

                sce.file = paste0(comp,".",status,".sce.RData")

                #check that SCE file exists (not missing due to missingness of ligand or receptor cell type)
                if(sce.file %in% basename(sce.files)){

                    #add ampersand or wait command, if batch size is met or last comparison
                    if( 
                       ((comp %% batch.size) == 0) 
                       | 
                       (comp == max.comps)
                    ){  
                        tail = "; wait;"
                    } else {
                        tail = " &"
                    }
          
                    #liana creates an "omnipath" temp folder using timestamps with seconds, so pause to avoid overwritting
                    cmd = paste0("sleep 2; Rscript ",comp,".",sample,".script.R &> ",comp,".",sample,".script.R.log",tail)

                    print(cmd)
                    cmds = c(cmds,cmd)
                }
            }
        }
    }
}

#save commands to file and run
cmds.file = paste0(comparisons.dir,"/liana.cmds.sh")

write.table(cmds, file = cmds.file, row.names=F,col.names=F,quote=F)

system(
    paste0(
        "sh ",cmds.file
    )
)

