#https://github.com/naoto-hikawa/RNAseq_tutorial/blob/main/DEG_tutorial/clean_to_DEG.R

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(dplyr)
library(tidyverse)


###########################################################################################
###DESeq2  DEseq2针对有生物学重复的样本。
###数据来源https://4va.github.io/biodatasci/r-rnaseq-airway.html#data_needed  
counts <- read.csv("airway_scaledcounts.csv")
metadata <-  read.csv("airway_metadata.csv")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)

dds <- DESeq(dds)

resultsNames(dds)

res1 <- results(dds,tidy=TRUE, name="dex_treated_vs_control")



data_geneid=read.csv("annotables_grch38.csv")
data_geneid=data_geneid[,c(1,3)]
# Create normalized count data file
DESeq2_normalize <- counts(dds, normalized=TRUE) 

#DESeq2_normalize=merge(data_geneid,DESeq2_normalize,by=1)
deg_normalize=merge(data_geneid,res1,by.x=1,by.y=1)


m_normalize=merge(data_geneid,DESeq2_normalize,by.x=1,by.y=0)

###res1 <- cbind(data_geneid, res1)
normalized_file=paste(name, "_normalized.txt", sep="")
write.table(DESeq2_normalize, file=paste("results/", normalized_file, sep=""), sep="\t", append=F, quote=F, row.names=F) 


################################################################################
metadata$dex=c("a","a","a","b","b","c","c","c")

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)

dds <- DESeq(dds)
resultsNames(dds) 
res1 <- results(dds, name="dex_c_vs_a") 
res1=as.data.frame(res1)
res2 <- results(dds, name="dex_b_vs_a") 
res2 <- results(dds, name="dex_b_vs_c") 


dds1 <- rlogTransformation(dds) 
exprSet_new=assay(dds1)



##=========================================================================================================
###edgeR  edgeR对于单个样本比较好
##BiocManager::install("edgeR")

library(edgeR)
colnames(counts)
count1=column_to_rownames(counts,"ensgene")
count2=count1[,c(1,3,5,7,2,4,6,8)]

group.list=c(rep("control",4),rep("treated",4))
group.list=factor(group.list)
group.list=relevel(group.list,ref = "control")

dge <- DGEList(counts=count2,group=group.list)

dge <- calcNormFactors(dge) 
design <- model.matrix(~0+group.list)
design1 <- model.matrix(~group.list)


rownames(design)<-colnames(dge)
colnames(design)<-levels(group.list)

dge <- estimateDisp(dge, design, robust = TRUE)

plotBCV(dge)


fit <- glmQLFit(dge, design, robust=TRUE)   #拟合模型
cntr.vs.T <- makeContrasts(control-treated, levels=design)
res <- glmQLFTest(fit, contrast=cntr.vs.T)   #统计检验


edger.deg=topTags(res, n = nrow(res$counts))
View(edger.deg[["table"]])


fit2 <- glmFit(dge, design)

fit2 <- glmLRT(fit2, contrast=c(-1,1)) 

DEG=topTags(fit2, n=nrow(exp))
DEG=as.data.frame(DEG)


########################################################

cpm <- cpm(count1)
lcpm <- cpm(count1, log=TRUE, prior.count=2)


##===========================================================================================
###  Linear Models for Microarray Data
library(limma)




group.list=c(rep("control",4),rep("treated",4))
group.list=factor(group.list)
group.list=relevel(group.list,ref = "control")


design <- model.matrix(~0+group.list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group.list)



dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")

fit <- lmFit(v, design)

constrasts = paste(rev(levels(group.list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
DEG = topTable(fit2, coef=constrasts,sort.by = "P", n=Inf)
DEG2=topTable(fit2, coef=1, adjust="BH")
DEG = na.omit(DEG)

  
  
  
  