#Summary: Generate heatmaps of donor assignment probabilities from Vireo
#Conda: cdh


library(pheatmap)

set.seed(100)
options(width=100)

source("utilities.r")


donors.dir = "../OutputData/Donors"


donor.order.map = list()
donor.order.map[["E215_NitrofenControl"]] = c(1,0,2,3,4)
donor.order.map[["E215_NitrofenCDH"]] = c(0,3,4,1,2)
donor.order.map[["New_E215_HealthyControl"]] = c(2,3,1,0,4)
donor.order.map[["Original_E215_HealthyControl"]] = c(2,3,1,0,4)


donors.dirs = list.dirs(path = donors.dir, full.names = T, recursive = F)

donors.res = foreach(dir = donors.dirs) %do% {
    print(dir)

    sample = basename(dir)

    vireo.path = paste0(dir,"/vireo.out")
    vireo.file = paste0(vireo.path,"/prob_singlet.tsv")
    cmd = paste0("zcat ",vireo.file,".gz > ",vireo.file)
    system(cmd)

    probs = read.table(vireo.file,header=T,sep="\t",row.names=1)

    pdf(file = paste0(vireo.path,"/", sample," - Vireo Donor Assignment Probabilities.pdf"))
    pheatmap(
        probs[,donor.order.map[[sample]]],
        show_rownames = F,
        cluster_cols = F,
        color = colorRampPalette(c("blue","yellow"))(100),
        scale = "none",
        fontsize_col = 20,
        main = ""
    )
    dev.off()
}

