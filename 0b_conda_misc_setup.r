#Summary: Additional code setup 

Sys.setenv(CC="x86_64-conda_cos6-linux-gnu-cc")
Sys.setenv(CXX="x86_64-conda_cos6-linux-gnu-c++")
Sys.setenv(CXX="x86_64-conda-linux-gnu-c++")

#libraries needed by Seurat
libs = read.csv("R_packages.txt",header=F,comment.char="#")$V1
libs

for(lib in libs){
    print(lib)
    
    if(lib %in% rownames(installed.packages()) == FALSE) {
        BiocManager::install(lib)
    }
    library(lib, character.only = T)
}


install.packages(
    pkgs = 'https://cran.r-project.org/src/contrib/spatstat.sparse_3.0-1.tar.gz', 
    repos=NULL
)


install.packages(
    pkgs = 'https://cran.r-project.org/src/contrib/Archive/spatstat.core/spatstat.core_2.4-4.tar.gz', 
    repos=NULL
)


#Seurat
install.packages(
    pkgs = "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.0.3.tar.gz",
    repos = NULL
)

