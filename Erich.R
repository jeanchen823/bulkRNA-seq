Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(dplyr)
library(tidyverse)

counts <- read.csv("airway_scaledcounts.csv")
metadata <-  read.csv("airway_metadata.csv")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)

dds <- DESeq(dds)
res1 <- results(dds,tidy=TRUE)

anno <- read.csv("annotables_grch38.csv")
anno=anno[,c(1,3)]
res2=merge(anno,res1,by=1)

res3=distinct(res2,symbol,.keep_all = T)

res3=na.omit(res3)

library(msigdbr)
library(dplyr)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(data.table)

log2(1.5)

select.FPKM <- (res3$baseMean > 100)
table(select.FPKM)
select.log2FC <- abs(res3$log2FoldChange) >1
table(select.log2FC)
select.qval <- (res3$padj < 0.05)
table(select.qval)

select.vec=select.FPKM & select.log2FC & select.qval 

table(select.vec)

degs.list=as.character(res3$symbol)[select.vec]
keytypes(org.Hs.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
dotplot(erich.go.BP,showCategory = 15)
barplot(erich.go.BP,showCategory = 8)
erich.go.BP=erich.go.BP@result
write.table(erich.go.BP,"erich.go.BP.deg.con.hf.txt",sep = "\t",col.names = NA)


erich.go.CC = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)


erich.go.MF = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.MF,showCategory = 8)

barplot(erich.go.MF,showCategory = 8)



####  

kegg=read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
length(unique(kegg$term))
res=enricher(gene = degs.list,TERM2GENE=kegg)

####################################################################################
##### GSEA
select.FPKM <- (res2$baseMean > 10)
table(select.FPKM)
rownames(res2)
degs.list=as.character(res2$symbol)[select.FPKM]
res3=distinct(res2,symbol,.keep_all = T)
rownames(res3)=res3$symbol
f1=res3[degs.list,]
exp.fc=f1[,c(2,4)]

### 

exp.fc.id.sorted <- exp.fc[order(exp.fc$log2FoldChange, decreasing = T),]
head(exp.fc.id.sorted)

id.fc.1 <- exp.fc.id.sorted$log2FoldChange
names(id.fc.1) <- exp.fc.id.sorted$symbol
kegg=read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
res= GSEA(id.fc.1,TERM2GENE=kegg)

gseaplot2(res, color = "red", pvalue_table = T,geneSetID = 1)

gseaplot2(res, geneSetID = c(2,3))


###################################################

deg = exp.fc.id.sorted$log2FoldChange
names(deg) = exp.fc.id.sorted$symbol
deg= sort(deg,decreasing = T)

library(msigdbr)
h.kegg= msigdbr(species ="Homo sapiens", category = "C2", subcategory = "KEGG")
length(unique(h.kegg$gs_exact_source))
h.kegg= h.kegg[,c(3,4)]
gsea <- GSEA(deg,pvalueCutoff = 0.5, TERM2GENE =h.kegg)
gseaplot2(gsea, geneSetID = 6, title = gsea$Description[6])
gsea1=gsea@result

write.table(gsea1,"gsea.res.txt",sep = "\t",col.names = NA)

#####################################################################################
a=degs.list[1:40]

b<-a[-which(is.na(a))]
write.table(a,file="deg.gene_symbol.txt",
            row.names = F,quote = F,col.names = F)
getwd()

##https://string-db.org/cgi/network?pollingId=bUpew6ZPM4MT&sessionId=bq3ONpHVnJ16

library(ggraph)
help(package="ggraph")
library(igraph)
nodes<-read.csv("deg.gene_symbol.txt",header=F)
links<-read.table("string_interactions_short(1).tsv",header=F,sep="\t")
nodes$Name<-nodes$V1
nodes
links
net<-graph_from_data_frame(d=links,vertices=nodes,directed = T)
plot(net)
ggraph(net,layout = "linear",circular=T)+
  geom_edge_link(color="blue")+
  geom_node_text(aes(label=Name))+
  geom_node_point(size=10,color="red",alpha=0.5)+
  theme_void()

?ggraph

## layout = circlepack
ggraph(net,layout = "treemap")+
  geom_edge_link(color="blue")+
  geom_node_text(aes(label=Name))+
  geom_node_point(size=10,color="red",alpha=0.5)+
  theme_void()

library(ggraph)
library(igraph)
nodes$Name<-nodes$V1
nodes$Group<-c(rep("Up",20),rep("Down",20))
net<-graph_from_data_frame(d=links,vertices=nodes,directed = T)
ggraph(net,layout = "kk")+
  geom_edge_link()+
  geom_node_point(size=10,aes(color=Group))+
  geom_node_text(aes(label=Name))+
  theme_void()

