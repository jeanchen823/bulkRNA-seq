Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())


library(GEOquery)
library(dplyr)
library(tidyverse)

#########################################################################################
gset = getGEO('GSE149507', destdir=".",getGPL =F)
soft=getGEO(filename ="GSE149507_family.soft")
gpl= soft@gpls[["GPL23270"]]@dataTable@table

gset=gset[[1]]
pdata=gset@phenoData@data

exprSet= gset@assayData[["exprs"]] 

#############################################

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
x1=gpl$ENTREZ_GENE_ID
x1=as.character(x1)
x2=AnnotationDbi::select(org.Hs.eg.db, keys=x1, 
                         columns=c("ENTREZID", "SYMBOL"), 
                         keytype="ENTREZID")
gpl=gpl[,c(1,2)]
gpl$gene=x2$SYMBOL


colnames(gpl)
exp=as.data.frame(exprSet)
exp.pl=merge(gpl,exp,by.x=1,by.y=0)
#############################################################################
rownames(exp.pl)=exp.pl$gene
colnames(exp.pl)
exp.pl.1=distinct(exp.pl,gene,.keep_all = T)
exp.pl.2=na.omit(exp.pl.1)
rownames(exp.pl.2)=exp.pl.2$gene
exp.pl.3=exp.pl.2[,-c(1:3)]

######################################################################################
a=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35)
a2=a+1
exp.4=exp.pl.3[,c(a2,a)]
p.n=filter(pdata,source_name_ch1=="normal")
p.t=filter(pdata,source_name_ch1=="tumor")
exp.4.1=exp.pl.3[,c(rownames(p.n),rownames(p.t))]

x=arrange(pdata,source_name_ch1)
exp.4.2=exp.pl.3[,rownames(x)]

rep('control',18)
group_list=c(rep('control',18),rep('T',18))
group_list=factor(group_list)
group_list <- relevel(group_list, ref="control")

library(limma) 
##除去批次效应
exp.5=normalizeBetweenArrays(exp.4)

boxplot(exp.4,outline=FALSE, notch=T,col=group_list, las=2)
boxplot(exp.5,outline=FALSE, notch=T,col=group_list, las=2)
library(RColorBrewer)
colors<-brewer.pal(18,"Dark2")
boxplot(exp.5,col=colors,notch=T,outline=FALSE, las=3,ylim=c(2,8))
boxplot(exp.4,col=colors,notch=T,outline=FALSE, las=3,ylim=c(2,10))


##构建差异分析的矩阵
design=model.matrix(~ group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp.5)


fit=lmFit(exp.5,design)
fit=eBayes(fit) 
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf) 
write.table(allDiff,file = "NA.allDiff.txt",sep = "\t",col.names = NA)
write.table(allDiff,file = "allDiff.txt",sep = "\t" )

write.csv(allDiff,file = "deg.csv")
###############################################################

###########################################################################################
##火山图

library(ggplot2)                         

library(ggrepel)                               

data <- allDiff
data$significant="stable"
log2(2)

data$significant[data$logFC>=1.5 & data$adj.P.Val <0.05]="up"

data$significant[ ]="up"

x=data$logFC>=1.5 & data$adj.P.Val <0.05
table(x)
x[1:20]

data$significant[data$logFC<= -1.5 & data$adj.P.Val <0.05]="down"

table(data$significant)

log10(0.05)

ggplot(      data,aes(       x=logFC,y=-1*log10(adj.P.Val)     )    )



ggplot(data,aes(x=logFC,y=-1*log10(adj.P.Val)))+xlim(-8,8)+ylim(0,12)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#f8961e","#EE7AE9"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.8)+
  geom_vline(xintercept = c(-1.5,1.5),linetype=4,size=0.8)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2 fold change ",y="-log10(p_value)")

select.FPKM <-  data$AveExpr > 5  
table(select.FPKM)

select.log2FC <- abs(data$logFC) >2
table(select.log2FC)
select.qval <-  data$adj.P.Val< 0.01 
table(select.qval)

select.vec= select.FPKM & select.log2FC & select.qval  
table(select.vec)

degs.list=as.character(rownames(data))[select.vec]

label.deg=sample(degs.list,20)

p=ggplot(data,aes(logFC,-1*log10(P.Value)))+xlim(-5,5)+ylim(0,10)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#EE7AE9","#f8961e"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-1.5,1.5),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")

p

x=arrange(data,logFC)

top.down=rownames(x)[1:10]
top.up=rownames(x)[20239:20251]
data_selected <- data[c(top.down,top.up),]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))

p + geom_text_repel(data=data_selected,
                    aes(label=rownames(data_selected)))
###########################

data_selected <- data[ c("ASCL1","SCG3","WIF1") ,]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))


####################################################3
##热图

annotation_col1 = data.frame(
  Database =c(rep("GEO",36)),
  CellType =c(rep("N",18),rep("T",18)) 
)
rownames(annotation_col1)=colnames(exp.5)
exp.6=as.data.frame(exp.5)
exprSet.map=exp.6[c(top.down,top.up),]
pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows = F,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   show_colnames=F,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("blue", "white","red"))(100))
####聚类
pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows =T,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   show_colnames=F,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("blue", "white","red"))(100))

######################################################
###箱线图
label.deg[1:10]
exp.x=exp.6["TTC32",]

z=t(exp.x)
z=as.data.frame(z)
rep("huage",10)

z$type=c(rep("N",18),rep("T",18))

colnames(z)
library(ggpubr)
library(ggsignif)
library(ggplot2)


ggboxplot(z,x="type",y="TTC32",
          width = 0.6,fill="type",
          notch = F,palette = c("#00AFBB", "red","#E7B800"),
          add = "jitter",shape="type")

###加p值
p=ggboxplot(z,x="type",y="TTC32",
            width = 0.6,fill="type",
            notch = F,palette = c("#00AFBB", "red","#E7B800"),
            add = "jitter",shape="type")
p
p + stat_compare_means(aes(group =type))

####################################################
###箱线图
ggplot(z,aes(type,TTC32,fill=type)) +
  geom_boxplot()+theme_classic()+
  xlab("huage")+ylab("exp")+ 
  geom_jitter(shape=9, position=position_jitter(0.2))

ggplot(z, aes(type,TTC32,fill=type)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5)

#######匹配
z$pairinfo=rep(1:18,2)
ggplot(z, aes(type,TTC32,fill=type)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of ","TTC32"))+
  theme_classic()+
  theme(legend.position = "none")

##########################################################################################
#####小提琴图

ggplot(z,aes(type,TTC32,fill=type)) +
  theme_classic()+
  geom_violin()+ geom_jitter(shape=16, position=position_jitter(0.2))

library(ggsignif)
compaired <- list(c("N", "T"))


library(ggsci)
colnames(z)
ggplot(z,aes(type,TTC32,fill=type)) +
  geom_violin()+
  theme_classic()+
  geom_signif(comparisons =compaired ,step_increase = 0.1,
              map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
              test =t.test,size=2,textsize = 6)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values = c("#2a9d8f", "#f8961e"))+geom_boxplot()

