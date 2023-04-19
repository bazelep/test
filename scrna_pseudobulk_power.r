#Summary: Calculate detectable fold-changes for endothelial cell genes based on pseudobulk counts.
#Conda: cdh



library(RNASeqPower)


#estimated number of genes with adequate expression to be included in DGE testing
n.genes = 10000 

#number of endothelial cell subtypes to test
n.cell.types = 3

#coefficient of variation for rat
cv = 0.1

#alpha adjusted for number of tests (#genes * #cell types)
alpha = 0.05 / (n.genes * n.cell.types)

#power
power = 0.8

#size of group 1, e.g. nitrofen control
n1 = 5

#size of group 2, e.g. nitrofen cdh
n2 = 5

#range of proportions of endothelial cells out of total cells
props = seq(1,50,by=2)/100

#total number of reads sequenced per sample
total = 200000000

#reads per donor
per.donor = total/n1

#sequencing depths (gene counts) for each donor for endothelial cells, adjusted for total number of genes assayed,
# for a range of EC proportions
depths = (props*per.donor)/25000

#calculate detectable fold-changes based on range of depths 
fc = rnapower(
    depth=depths,
    cv=cv,
    alpha=alpha,
    power=power,
    n=n1,
    n2=n2
)
              

results = data.frame("EC.Proportion"=props,"Detectable.Fold.Change"=fc)

write.csv(results,file="scRNA.EndothelialCell.Donor.Pseudobulk.DetectableFoldChanges.csv")

