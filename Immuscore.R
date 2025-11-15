Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(GEOquery)
library(dplyr)
library(tidyverse)
library(data.table)

gset = fread("GSE19188_series_matrix.txt",data.table = F,sep = "\t",skip=68)

anno=fread("GPL570-55999.txt",sep = "\t",header = T,data.table = F)

gene=anno[,c(1,11)]


exp1=as.data.frame(gset)
exp.anno=merge(x=gene,y=exp1,by.x=1,by.y=1)
x=gene$`Gene Symbol`

a1=strsplit(x,split = " /// ",fixed = T)

gene.all = sapply(a1,function(x){x[1]})
gene$gene.all=gene.all
gene=distinct(gene,gene.all,.keep_all = T)

gene1=na.omit(gene)

rownames(exp.anno)=exp.anno$ID
exp.anno1=exp.anno[gene1$ID,]
exp.anno2=exp.anno1[,-c(1,2)]
rownames(exp.anno2)=gene1$gene.all

write.table(exp.anno2,"exp.anno2.txt",sep = "\t",col.names = NA)
##BiocManager::install("tidybulk")
source("Cibersort.R")
ciber <- CIBERSORT(sig_matrix = "LM22.txt",
                   mixture_file = "exp.anno2.txt",
                   perm = 100,
                   QN = TRUE)

phe=fread("phe.txt",data.table = F)
colnames(phe)[1]="name"
phe$name= paste0("name", 1:37)
rownames(phe)=phe$name
phe=phe[,-1]
phe=t(phe)
phe=as.data.frame(phe)
table(phe$name8)
t.luad=filter(phe,name8=="tissue type: tumor")
ciber1=as.data.frame(ciber)
ciber.t=ciber1[rownames(t.luad),c(1:22)]
exp.anno3=exp.anno2[,rownames(t.luad)]

tcga_expr=exp.anno3
tcga_gsva=ciber.t

immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]),y,type="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

CD274=immuscore("CD274")
IFNG=immuscore("IFNG")
CD274$dataset="GSE19188"

CD274$zf=ifelse(CD274$cor>0,"Positive","Negative")
CD274$cor1=abs(CD274$cor)
save(CD274,file = "CD274.g19188.rdata")


colnames(CD274)
ggplot(CD274,aes(x=dataset,y=immune_cells))+
  geom_point(aes(color=zf,size=cor1),alpha=0.8)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5))
