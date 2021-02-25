if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(limma)
library(affy)
library(hgu133acdf)


#IRF4 - is there a difference in gene expression between mBl and intermediate


setwd("C:/Users/Kai/Desktop/IRF4/Raw")
batch = read.affybatch(dir(patt="CEL")) # puts all .cel file info into one thing

IRF4sdrf = read.csv(file = "C:/Users/Kai/Desktop/IRF4/IRF4sdrf.csv", header=T)

norm.batch = rma(batch) # we perform RMA normalisation on our expression data across all samples.
pData(norm.batch) # data frame which shows the sample's name and number in norm batch
dat = exprs(norm.batch)
print(dat)

t2 = vector()
pval.t2 = vector()
group = IRF4sdrf$Characteristics.molecular.diagnosis.
for(j in 1:nrow(dat)){
  temp = dat[j,]
  res=t.test(temp[group=="intermediate"], temp[group=="non-mBL"], var.equal=T)
  t2[j] = res$stat
  pval.t2[j]=res$p.val
}



adj.pval.t2 <- p.adjust(pval.t2, "BH")
result.table2 = data.frame(ID=rownames(dat), t.stat=t2,
                           pvalue=pval.t2, fdr.pvalue=adj.pval.t2)
result.table2.sorted = result.table2[order(adj.pval.t2),]
results <-result.table2.sorted[1:12,] # listing the top 10 genes
results

library("annotate")
BiocManager::install("hgu133a.db")
library("hgu133a.db")
GENEID <-select(hgu133a.db, c("203434_s_at", "203435_s_at", "211862_x_at", "215731_s_at","219515_at",
                     "203680_at", "202716_at","203244_at","204490_s_at","212014_x_at", "213899_at", "222270_at"), c("SYMBOL","ENTREZID", "GENENAME")) #getting the gene id's for the top 10
GENEID

GENEIDs <- select(hgu133a.db, c(result.table2.sorted$ID), c("SYMBOL","ENTREZID", "GENENAME")) # gene Id's for all genes
GENEIDs

GeneIDs <- as.data.frame(GENEIDs)
write_xlsx(GeneIDs, "C:\\Users\\Kai\\Desktop\\IRF4\\samrIDs.xlsx")


#SAM 

BiocManager::install("samr")
install.packages("samr", repos= "	https://statweb.stanford.edu/~tibs/SAM", type = "source")
library(samr)

install.packages("writexl")
library(writexl)

write_xlsx(samrIDs, "C:\\Users\\Kai\\Desktop\\IRF4\\samrIDs.xlsx") # probe IDs to gene symbols, merged to data 
dat <- as.data.frame(dat)
write_xlsx(dat, "C:\\Users\\Kai\\Desktop\\IRF4\\RMAsam.xlsx") #data itself

#opening a seperate window to conduct SAM analysis

runSAM()

SAMID <-select(hgu133a.db, c("212014_x_at","209835_x_at","212063_at","211316_x_at"), c("SYMBOL","ENTREZID", "GENENAME"))
#getting the gene id's for the top 10

SAMID


















