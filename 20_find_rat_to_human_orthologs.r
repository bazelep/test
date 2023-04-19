#Summary: Find rat to human orthologous genes to use with ligand-receptor analysis.
#Conda: cdh


library(biomaRt)

set.seed(100)
options(width=100)



counts.dir = "../OutputData/Counts"
ligand.receptor.dir = "../OutputData/Receptor.Ligand"

if(!dir.exists(ligand.receptor.dir)){
    dir.create(ligand.receptor.dir)
}


#get full set of ensembl ids since more reliable than gene symbols
# - pick any sample
cmd = paste0("zcat ",counts.dir,"/E21_5_LeftLungControl/outs/filtered_feature_bc_matrix/features.tsv.gz | awk -F\"\t\" '{print $1}'")
ensembl.gene.ids = system(cmd,intern=T)


#create biomart connection
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <- useEnsembl(biomart = "ensembl", version = 108)
ensembl.rat <- useDataset(dataset = "rnorvegicus_gene_ensembl", mart = ensembl)

#query for human orthologs
features.orthologs = getBM(
      attributes= c("ensembl_gene_id",
                    "hsapiens_homolog_ensembl_gene",
                    "hsapiens_homolog_associated_gene_name",
                    "external_gene_name"),
      filters=c("ensembl_gene_id"),
      values=ensembl.gene.ids,
      mart=ensembl.rat)

names(features.orthologs) = c("Rat.Ensembl.Gene.ID","Human.Ensembl.Gene.ID","Human.gene.name","Rat.gene.name")

saveRDS(features.orthologs, file = paste0(ligand.receptor.dir,"/rat_to_human_orthologs.rds"))

