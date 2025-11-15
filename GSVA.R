Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(data.table)
library(dplyr)
library(tidyverse)

stad.phe=fread("TCGA-STAD.GDC_phenotype.tsv.gz",header = T, sep = '\t',data.table = F)

stad.fkpm=fread("TCGA-STAD.htseq_fpkm.tsv.gz",header = T, sep = '\t',data.table = F)

stad.pro=fread("gencode.v22.annotation.gene.probeMap",
               header = T, sep = '\t',data.table = F)
stad.pro=stad.pro[,c(1,2)]

stad.fkpm.pro=merge(stad.pro,stad.fkpm,by.y ="Ensembl_ID",by.x = "id" )
dim(stad.fkpm.pro)
distinct()
stad.fkpm.pro[,-c(1,2)]
stad.fkpm.pro.1 <- aggregate(stad.fkpm.pro[,-c(1,2)], list(stad.fkpm.pro$gene), FUN=sum)

####################################################################
stad.fkpm.pro=distinct(stad.fkpm.pro,gene,.keep_all = T)

stad.fkpm.pro=distinct(stad.fkpm.pro,gene,.keep_all = T)
dim(stad.fkpm.pro)
#rownames(stad.fkpm.pro)=stad.fkpm.pro$gene
#stad.fkpm.pro=stad.fkpm.pro[,-2]
stad.fkpm.pro <- column_to_rownames(stad.fkpm.pro,"gene")
rownames(stad.phe)=stad.phe$submitter_id.samples
colnames(stad.phe)
table(stad.phe$sample_type.samples)
stad.phe.t=filter(stad.phe,sample_type.samples=="Primary Tumor")
stad.phe.n=filter(stad.phe,sample_type.samples=="Solid Tissue Normal")
intersect(c("a","b"),c("b","c"))
z1=intersect(rownames(stad.phe.t),colnames(stad.fkpm.pro))
z2=intersect(rownames(stad.phe.n),colnames(stad.fkpm.pro))
stad.t=stad.fkpm.pro[,z1]
stad.n=stad.fkpm.pro[,z2]
colnames(stad.n)=paste0("N",1:32)
paste("xy",c(1,5,8,11,3),sep = "-")
?paste
colnames(stad.t)=paste0("T",1:375)
stad.exp=merge(stad.n,stad.t,by.x = 0,by.y = 0)
colnames(stad.exp)
stad.exp <- column_to_rownames(stad.exp,"Row.names")
gtex.exp=fread("gtex_RSEM_gene_fpkm.gz",header = T, sep = '\t',data.table = F)

gtex.exp[1:5,1:5]

gtex.phe=fread("GTEX_phenotype.gz",header = T, sep = '\t',data.table = F)

gtex.pro=fread("probeMap_gencode.v23.annotation.gene.probemap",header = T, sep = '\t',data.table = F)
gtex.pro=gtex.pro[,c(1,2)]
###gtex.exp=merge(gtex.exp,gtex.pro,by.x ="sample",by.y = "id" )

######################################################################
save(stad.exp,file = "stad.exp.rdata")
rm(stad.exp,stad.fkpm,stad.fkpm.pro,stad.n,stad.t,stad.phe)


colnames(gtex.phe)
rownames(gtex.phe)=gtex.phe$Sample
table(gtex.phe$primary_site)
colnames(gtex.phe)[3]= "primary_site" 
colnames(gtex.phe)
table(gtex.phe$primary_site)
gtex.phe.s=filter(gtex.phe,primary_site=="Stomach")
x1=intersect(rownames(gtex.phe.s),colnames(gtex.exp))
colnames(gtex.exp)[1:10]
gtex.s=gtex.exp[,c("sample",x1)]

gtex.exp=merge(gtex.pro,gtex.s,by.x ="id",by.y ="sample")


gtex.s1=distinct(gtex.exp,gene,.keep_all = T)
gtex.s2 <- column_to_rownames(gtex.s1,"gene")
gtex.s2 =gtex.s2[,-1]


gtex.s3=2^gtex.s2

gtex.s5=log2(gtex.s3-0.001+1)
colnames(gtex.s5)=paste0("G",1:174)

length(intersect(gtex.exp$id,stad.fkpm.pro$id))
length(intersect(gtex.exp$gene,rownames(stad.fkpm.pro)))

all.data=merge(gtex.s5,stad.exp,by.x = 0,by.y = 0)
all.data <- column_to_rownames(all.data,"Row.names")
fwrite(all.data,file = "stomach.cancer.gtex.tcga.txt",sep = "\t",row.names = T)
write.table(all.data,file = "stomach.cancer.gtex.tcga.txt",sep = "\t")
all.data1<- as.matrix(all.data)

#################################################################################################

options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("limma")
library(limma)
nromalized.data=normalizeBetweenArrays(all.data)


library(GSVA)
library(GSEABase)
##BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install.packages("devtools")
#devtools::install_github("GSEA-MSigDB/GSEA_R")
#library(GSEA)

kegggmt2 <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")

kegg_list = split(kegggmt2$gene, kegggmt2$term)

##BiocManager::install("GSVA")
library(GSVA)

nromalized.data=as.matrix(nromalized.data)
kegg2 <- gsva(nromalized.data, kegg_list, kcdf="Gaussian",method = "gsva")



annotation_col = data.frame(
  Tissuetype =c(rep("Stomach",174),rep("Solid Tissue Normal",32),rep("Tumor",375)),
  Database =c(rep("GTEX",174), rep("TCGA",407))
)

rownames(annotation_col)=colnames(nromalized.data)

pheatmap::pheatmap(kegg2[20:40,], #热图的数据
                   cluster_rows =T,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   annotation_col = annotation_col,
                   show_colnames=F,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))


###############################################################################################
##不同物种

##install.packages("msigdbr")
library(msigdbr)
msigdbr_species()

h.kegg= msigdbr(species ="Homo sapiens", category = "C2", subcategory = "KEGG")

h.kegg=h.kegg[,c(3,4)]
h.kegg_list = split(h.kegg$gene_symbol, h.kegg$gs_name)

###############################################################
#### https://github.com/Horizondmh/R-script-for-scRNA-seq-analysis/blob/main/scNMO_CART.R


expr=nromalized.data
human_GO_bp_Set=kegg_list 
gsva.res = gsva(expr, human_GO_bp_Set, method = "ssgsea")
gsva_gobp = data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva_gobp, "B_gsva_GOBP.csv", row.names = F)
library(limma)
gsva_gobp = read.csv("B_gsva_GOBP.csv",row.names = 1)

#Tissuetype =c(rep("Stomach",174),rep("Solid Tissue Normal",32),rep("Tumor",375)),
#Database =c(rep("GTEX",174), rep("TCGA",407))

group <- factor(c(rep("Normal",206), rep("Tumor", 375)), levels = c('Normal', 'Tumor'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(gsva_gobp)
design
# PBMC VS CSF
compare <- makeContrasts(Normal - Tumor, levels=design)
fit <- lmFit(gsva_gobp, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1,n=Inf)
head(Diff)
write.csv(Diff,"diff_B_GOBP.csv")
library(stringr)
Diff$id = rownames(Diff)
Diff$id <- str_replace(Diff$id , "KEGG_","")
Diff$id <- gsub("_"," ",Diff$id)
Diff$id <- tolower(Diff$id)
write.csv(Diff,"diff_PPB.csv")
Diff =  Diff[c(1:6,41:45),]
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
library(ggplot2)
library(ggthemes)
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= "#D24B27",'NoSignifi'='#cccccc','Down'="#3BBCA8")) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, B_PBMC versus B_CSF') +
  guides(fill=F)+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
low1 <- dat_plot %>% filter(t < -10) %>% nrow()
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 10) %>% nrow()
high1 <- nrow(dat_plot)
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') +
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + 
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + 
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 
p

########################################################################################

##自定义基因集、另外一种gmt


gene.set=read.table("GENE.txt",
                    header =F,sep = '\t',quote = '')
kegg.name=read.table("GSET.NAME.txt",
                     header =F,sep = '\t',quote = '')
gene.set1=as.matrix(gene.set)
gene.set2=t(gene.set1)
gmt=list()
for (i in 1:24) {
  y=as.character(gene.set2[,i]) 
  b=which(y=="")
  gmt[[i]]=y[-b]
}
names(gmt)=kegg.name$V1
gmt=gmt[-24]


###################################################

y=as.character(gene.set2[,1]) 
b=which(y=="")
gmt[[1]]=y[-b]

View(gmt)
getwd()





###############特定基因分析
exprSet.all.r=all.data[c("RORC", "RORB", "RORA"),]

exprSet.all.r=t(exprSet.all.r)
exprSet.all.r=as.data.frame(exprSet.all.r)

x=c(rep("GTEX",174),rep("N",32),rep("T",375))
exprSet.all.r$Type=x

exprSet.rorc=exprSet.all.r[,c(1,4)]
exprSet.rorc$Gene=rep("RORC")
colnames(exprSet.rorc)[1]="Relative Expression"
exprSet.rorb=exprSet.all.r[,c(2,4)]
exprSet.rorb$Gene=rep("RORB") 
colnames(exprSet.rorb)[1]="Relative Expression"
exprSet.rora=exprSet.all.r[,c(3,4)] 
exprSet.rora$Gene=rep("RORA") 
colnames(exprSet.rora)[1]="Relative Expression"
x.all=rbind(exprSet.rorc,exprSet.rorb,exprSet.rora)
colnames(x.all)
library(ggsignif)
library(ggpubr)
library(ggplot2)
p <- ggboxplot(x.all, x = "Gene", y = "Relative Expression",
               color = "Type", palette = "Type",
               add = "Type")
p + stat_compare_means(aes(group = Type))

table(x.all$Gene)
my_comparisons <- list(c("RORA","RORB"), c("RORA","RORC"),c("RORB", "RORC"))


p +geom_signif(comparisons = my_comparisons,
               step_increase = 0.2,map_signif_level = F,
               test = t.test,size=0.8,textsize =4)
?geom_signif
x.c.b=cbind(exprSet.rorc,exprSet.rorb)
GSEA
GSVA
GO/KEGG
colnames(x.c.b)=c("RORC","Type","Gene", "RORB","Type","Gene" )
x.c.b=x.c.b[,c(1,4)]
library(ggplot2)
library(ggpubr)
## Loading required package: magrittr
p1 <- ggplot(data = x.c.b, mapping = aes(x = RORC, y = RORB)) +
  geom_point(colour = "red", size = 2) +
  geom_smooth(method = lm, colour='blue', fill='gray') #添加拟合曲线
p1
p1 + stat_cor(method = "pearson", label.x =3, label.y = 1) #添加pearson相关系数


##############################################################################################
#####生存曲线
stad.fkpm.pro.1=column_to_rownames(stad.fkpm.pro.1,colnames(stad.fkpm.pro.1)[1])

library(limma)
stad.fkpm.pro.1=normalizeBetweenArrays(stad.fkpm.pro.1)
kegg_list=kegg_list[1:10]
kegg2 <- gsva(stad.fkpm.pro.1, kegg_list,method = "gsva")


sur=fread("TCGA-STAD.survival.tsv",data.table = F)

sur1 <- sur[sur$sample %in% colnames(stad.fkpm.pro.1),]
##将TCGA表达矩阵的病人样本与临床数据样本进行整合

sur.kegg=as.data.frame(kegg2)

sur.kegg.1=sur.kegg[2,]
sur.kegg.2=as.data.frame(t(sur.kegg.1))
colnames(sur.kegg.2)="Glycan"

sur.glycan <- merge(sur1, sur.kegg.2, by.x="sample", by.y=0)


##中位数分组
sur.glycan$m = ifelse(sur.glycan$Glycan>median(sur.glycan$Glycan),'high','low')


##################生存分析

library(survival)
#install.packages("survminer")
library(survminer)


################################################################################
colnames(sur.glycan)
ggsurvplot(survfit(Surv(OS.time,OS)~m, data=sur.glycan), conf.int=T, pval=TRUE)

ggsurvplot(survfit(Surv(OS.time,OS)~m, data=sur.glycan), conf.int=F, pval=T)


fit=survfit(Surv(OS.time,OS)~m, data=sur.glycan)
ggsurvplot(fit,   
           title = "Survival curve", # 添加标题
           font.main = c(16, "bold", "darkblue"), # 设置标题字体大小、格式和颜色
           font.x = c(14, "bold.italic", "red"), # 设置x轴字体大小、格式和颜色
           font.y = c(14, "bold.italic", "darkred"), # 设置y轴字体大小、格式和颜色
           font.tickslab = c(12, "plain", "darkgreen")) # 设置坐标轴刻度字体大小、格式和颜色


############################################################################


###四分位生存曲线
quantile(sur.glycan$Glycan)
s2.exp$zz=ifelse(s2.exp$EZH2>9.53,'high', 
                 ifelse( s2.exp$EZH2>9.23& s2.exp$EZH2<9.53,'h.stable', 
                         ifelse( s2.exp$EZH2>8.88& s2.exp$EZH2<9.23,'l.stable','down') )
)


##########################################################################
