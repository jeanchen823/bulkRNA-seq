Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

##BiocManager::install("ALL")
library(ALL)

data(ALL)
d = exprs(ALL)

mads = apply(d, 1, mad)
d = d[rev(order(mads))[1:5000], ]
d = sweep(d, 1, apply(d, 1, median, na.rm = TRUE))

## BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(d, maxK = 5, reps = 30, pItem = 0.8, pFeature = 1,
                               title = title, clusterAlg = "hc", distance = "pearson", 
                               seed = 123456, plot = "png")
# 计算聚类共识和项目共识
icl = calcICL(results, title = title, plot = "png")

icl[["clusterConsensus"]]

icl[["itemConsensus"]][1:5,]

Kvec = 2:20
x1 = 0.1
x2 = 0.9

PAC = rep(NA, length(Kvec))
names(PAC) = paste("K=", Kvec, sep="")

results = ConsensusClusterPlus(d, maxK =20, reps = 30, pItem = 0.8, pFeature = 1,
                               title = title, clusterAlg = "hc", distance = "pearson", 
                               seed = 123456, plot = "png")

for (i in Kvec) {
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])  # 计算共识矩阵的经验累积分布函数
  PAC[i-1] = Fn(x2) - Fn(x1)
}

optK = Kvec[which.min(PAC)]

##BiocManager::install("ConsensusClusterPlus")
#install.packages("dendsort")
load("stad.exp.rdata")
library(ConsensusClusterPlus)
library(pheatmap)
library(dendsort)
library(dplyr)
expr <- stad.exp[rowSums(stad.exp) > 1,]
colnames(expr)
expr=expr[,33:407]
library(msigdbr)

homo= msigdbr(species = "Homo sapiens")
REACTOME_GLYCOLYSIS.gene=filter(homo,gs_name=="REACTOME_GLYCOLYSIS")
REACTOME_CHOLESTEROL_BIOSYNTHESIS.gene=filter(homo,gs_name=="REACTOME_CHOLESTEROL_BIOSYNTHESIS")


purple <- "#3A236B"
lightblue <- "#A7DCF8"
seagreen <- "#7B9C4A"
nake <- "#BDA597"
cherry <- "#7F1D47"

glycolytic.gene <- REACTOME_GLYCOLYSIS.gene$gene_symbol # 基因数目与原文不符 
cholesterogenic.gene <- REACTOME_CHOLESTEROL_BIOSYNTHESIS.gene$gene_symbol

glycolytic.gene <- intersect(glycolytic.gene,rownames(expr))
cholesterogenic.gene <- intersect(cholesterogenic.gene,rownames(expr))

indata <- expr[c(glycolytic.gene,cholesterogenic.gene),] 
indata <- t(scale(indata))   # 数据标准化作聚类


cc <- ConsensusClusterPlus(d = indata, 
                           maxK = 10, # 三类
                           reps = 100, # 重抽样次数，可修改
                           pItem = 0.8, # 列维度抽样概率，可修改
                           pFeature = 1, # 行维度抽样概率，可修改
                           clusterAlg = "hc", # 聚类算法，可修改
                           innerLinkage = "ward.D2", # 内部链接函数，可修改
                           finalLinkage = "ward.D2", # 最终链接函数，可修改
                           distance = "spearman", # 距离测度，可修改
                           seed = 20240716,
                           title = "ConsensusCluster",
                           plot = "pdf") #或png


annCol <- data.frame(Cluster = paste0("C",cc[[3]][[3]]),
                     row.names = colnames(indata))
annRow <- data.frame(Pathway = rep(c("Glycolytic","Cholesterogenic"),c(length(glycolytic.gene),length(cholesterogenic.gene))),
                     row.names = colnames(indata))
annColors <- list(Cluster = c("C1" = seagreen,"C2" = cherry, "C3" = nake),
                  Pathway = c("Glycolytic" = "black","Cholesterogenic" = lightblue))

plotdata <- cc[[3]][[1]]
dimnames(plotdata) <- list(colnames(indata),colnames(indata))
pheatmap(mat = plotdata,
         color = colorRampPalette((c("white",purple)))(64),
         border_color = NA,
         cluster_rows = dendsort(cc[[3]][[2]]),
         cluster_cols = dendsort(cc[[3]][[2]]),
         annotation_col = annCol,
         annotation_row = annRow,
         annotation_colors = annColors,
         show_colnames = F,
         show_rownames = F)

dev.copy2pdf(file = "consensus heatmap.pdf",width = 6, height = 5)
?dev.copy2pdf


#绘制四象限散点图 - Figure 1B
glycolytic.curated <- setdiff(rownames(annCol[which(annCol$Cluster == "C2"),,drop = F]),cholesterogenic.gene) # 这里是C3
cholesterigenic.curated <- setdiff(rownames(annCol[which(annCol$Cluster == "C3"),,drop = F]),glycolytic.gene) # 这里是C2

pathway.list <- list("glycolytic" = glycolytic.curated,
                     "cholesterigenic" = cholesterigenic.curated)

curated.meta.expr <-  expr[c(glycolytic.curated,cholesterigenic.curated),]  
curated.meta.expr <- as.data.frame(t(scale(t(curated.meta.expr))))
curated.meta.expr$Pathway <- rep(c("glycolytic","cholesterigenic"),c(length(glycolytic.curated),length(cholesterigenic.curated)))
curated.meta.expr.median <- apply(curated.meta.expr[,setdiff(colnames(curated.meta.expr), "Pathway")], 2, function(x) tapply(x, INDEX=factor(curated.meta.expr$Pathway), FUN=median, na.rm=TRUE))
curated.meta.expr.median <- as.data.frame(t(curated.meta.expr.median))

curated.meta.expr.median$subtype <- 
  ifelse(curated.meta.expr.median$cholesterigenic >= 0 & curated.meta.expr.median$glycolytic >= 0,"Mixed",
         ifelse(curated.meta.expr.median$cholesterigenic <= 0 & curated.meta.expr.median$glycolytic <= 0,"Quiescent",
                ifelse(curated.meta.expr.median$cholesterigenic > 0 & curated.meta.expr.median$glycolytic < 0,"Cholesterogenic","Glycolytic")))

par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,4.1,2.1,3.1)+.1, las=1, tcl=-.25)
plot(NULL, NULL, # 绘制空白背景
     xlim = range(curated.meta.expr.median$glycolytic),
     ylim = range(curated.meta.expr.median$cholesterigenic),
     xlab = "median glycolytic gene expression (z-score)",
     ylab = "median cholesterigenic\ngene expression (z-score)")
grid(col = "grey85", lty = 2, lwd = 1.5) # 添加网格线
abline(v = 0, lty = 2, lwd = 2) # 添加水平0截断
abline(h = 0, lty = 2, lwd = 2) # 添加垂直0阶段

points(curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Mixed"),"glycolytic"],
       curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Mixed"),"cholesterigenic"],
       pch = 19,
       col = "#895430")

points(curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Cholesterogenic"),"glycolytic"],
       curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Cholesterogenic"),"cholesterigenic"],
       pch = 19,
       col = "#99B15B")

points(curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Glycolytic"),"glycolytic"],
       curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Glycolytic"),"cholesterigenic"],
       pch = 19,
       col = "#62A2AA")

points(curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Quiescent"),"glycolytic"],
       curated.meta.expr.median[which(curated.meta.expr.median$subtype == "Quiescent"),"cholesterigenic"],
       pch = 19,
       col = "#E6B73C")

legend("topleft",
       legend = c("Metabolic subgroup","Mixed","Cholesterogenic","Glycolytic","Quiescent"),
       col = c(NA,"#895430","#99B15B","#62A2AA","#E6B73C"),
       pch = 19,
       border = NA,
       bty = "n")

dev.copy2pdf(file = "median classification.pdf",width = 5, height = 5)

##绘制基于四分类的表达谱热图 - Figure 1C
plotdata <- curated.meta.expr[,setdiff(colnames(curated.meta.expr),"Pathway")] # 提出表达值
annCol <- data.frame("Subgroup" = curated.meta.expr.median$subtype,
                     row.names = rownames(curated.meta.expr.median),
                     stringsAsFactors = F)
annCol$Subgroup <- factor(annCol$Subgroup,levels = c("Quiescent","Glycolytic","Cholesterogenic","Mixed"))
annCol <- annCol[order(annCol$Subgroup),,drop = F]
annRow <- data.frame(Pathway = rep(c("Glycolytic","Cholesterogenic"),c(length(glycolytic.curated),length(cholesterigenic.curated))),
                     row.names = c(glycolytic.curated,cholesterigenic.curated),
                     stringsAsFactors = F)
annColors <- list(Subgroup = c("Mixed" = "#895430", "Cholesterogenic" = "#99B15B", "Glycolytic" = "#62A2AA", "Quiescent" = "#E6B73C"),
                  Pathway = c("Glycolytic" = "black","Cholesterogenic" = lightblue))

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}



plotdata <- standarize.fun(plotdata[rownames(annRow),rownames(annCol)],halfwidth = 4)
pheatmap(plotdata,
         border_color = NA,
         cluster_cols = F,
         cluster_rows = F,
         show_colnames = F,
         show_rownames = T,
         annotation_row = annRow,
         annotation_col = annCol,
         annotation_colors = annColors,
         gaps_row = 28, # 热图基因截断位置
         gaps_col =  cumsum(table(annCol$Subgroup))[1:3], # 热图亚型截断位置
         color = colorRampPalette((c("#18469B","white","#8C183D")))(64)) # 例文颜色




