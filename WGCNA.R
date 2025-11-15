library(data.table)
## BiocManager::install("WGCNA")
library(WGCNA)

library(tidyverse)
library(reshape2)

setwd("F:/shangke/lession11/s11/")
## BiocManager::install("GEOquery")
library(GEOquery)
geo=getGEO(filename='GSE48213_family.soft')
metadata=geo@gsms


meta1 = sapply(metadata,function(x){x@header[[2]]})
meta2=as.data.frame(t(meta1))



datTraits = data.frame(gsm=rownames(meta2),
                       cellline=trimws(sapply(as.character(meta2$V1),function(x) strsplit(x,":")[[1]][2])),
                       subtype=trimws(sapply(as.character(meta2$V3),function(x) strsplit(x,":")[[1]][2]))
)


getwd()
##########################################################################################
a = list.files()
n = length(a)

merge_data = read.csv(file = a[1],sep= "\t",header=T)
colnames(merge_data)
library(dplyr)
for (i in 2:n){
  
  new_data = read.csv(file = a[i],sep= "\t",header=T)
  
  merge_data = full_join(merge_data,new_data,by = "EnsEMBL_Gene_ID")
  
}

merge_data=column_to_rownames(merge_data,"EnsEMBL_Gene_ID")

#######################################################################################

merge_data = t(merge_data[order(apply(merge_data, 1, mad), 
                                decreasing = T)[1:5000],])

merge_data[1:4,1:4]

gene=merge_data

gsg = goodSamplesGenes(gene, verbose = 3)
gsg$allOK
gsg$goodSamples

###################################################################################################

powers <- 1:20

sft <- pickSoftThreshold(gene, powerVector = powers, verbose = 5)

par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type = 'n', 
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit,signed R^2', 
     main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, col = 'red');
abline(h = 0.90, col = 'red')

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n',
     main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')


sft$powerEstimate


###############################################################################################

power = sft$powerEstimate
power

#----------------------------------------------------------------------------------
#  Step 5: One-step network construction and module detection
#----------------------------------------------------------------------------------
net = blockwiseModules(gene, power = power, maxBlockSize = ncol(gene),
                       TOMType = "signed", minModuleSize = 20, 
                       reassignThreshold = 0, mergeCutHeight = 0.0001,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = "pearson", 
                       maxPOutliers = 1, loadTOMs=TRUE,
                       randomSeed = 931, # seed
                       saveTOMFileBase = paste0("F:/shangke/lession11/s11/"),
                       verbose = 3, pearsonFallback = "none", deepSplit = 3 )
# module:
table(net$colors) # 0 corresponds to unclassified genes
# Convert [number] labels to [colors] for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# plot
pdf(paste0( "./3.module.pdf"), width = 8, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


###############################################################################################
###############################################################################################

powers <- sft$powerEstimate

adjacency <- adjacency(gene, power = powers)
tom_sim <- TOMsimilarity(adjacency)
?TOMsimilarity
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)
tom_sim[1:6,1:6]

#write.table(tom_sim, 'TOMsimilarity.txt', sep = '\t', col.names = NA, quote = FALSE)

#####################################################################################################
#共表达模块识别

tom_dis  <- 1 - tom_sim


geneTree <- hclust(as.dist(tom_dis), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = FALSE, hang = 0.04)

minModuleSize <- 30  #模块基因数目
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05,
                    main = 'Gene dendrogram and module colors')



#基因表达聚类树和共表达拓扑热图
plot_sim <- -(1-tom_sim)
#plot_sim <- log(tom_sim)
diag(plot_sim) <- NA
TOMplot(plot_sim, geneTree, dynamicColors,
        main = 'Network heatmap plot, selected genes')


##########################################################################################


#计算基因表达矩阵中模块的特征基因（第一主成分）
MEList <- moduleEigengenes(gene, colors = dynamicColors)
MEs <- MEList$eigengenes
head(MEs)[1:6]
#write.table(MEs, 'moduleEigengenes.txt', sep = '\t', col.names = NA, quote = FALSE)

############################################################################################


#共表达模块的进一步聚类

ME_cor <- cor(MEs)
ME_cor[1:6,1:6]

#绘制聚类树观察
METree <- hclust(as.dist(1-ME_cor), method = 'average')
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')

abline(h = 0.4, col = 'blue')
abline(h = 0.5, col = 'red')

#模块特征基因聚类树热图
plotEigengeneNetworks(MEs, '', cex.lab = 0.8, xLabelsAngle= 90,
                      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2))

merge_module <- mergeCloseModules(gene, dynamicColors, cutHeight = 0.5, verbose = 3)
mergedColors <- merge_module$colors
table(mergedColors)

#基因表达和模块聚类树
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c('Dynamic Tree Cut', 'Merged dynamic'),
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05)


#################################################################################################################################################


#共表达模块与与临床表型的关联分析

module <- merge_module$newMEs
trait=datTraits

nGenes = ncol(gene)
nSamples = nrow(gene)
x=as.factor(datTraits$subtype) 
design=model.matrix(~0+ datTraits$subtype)
colnames(design)=levels(x)


trait=design

module <- merge_module$newMEs

moduleTraitCor <- cor(module, trait, use = 'p')
moduleTraitCor[1:6,1:5]  #相关矩阵

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(module))

#write.table(moduleTraitCor, 'moduleTraitCor.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.table(moduleTraitPvalue, 'moduleTraitPvalue.txt', sep = '\t', col.names = NA, quote = FALSE)

#相关图绘制
textMatrix <- paste(signif(moduleTraitCor, 2), '\n(', signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, main = paste('Module-trait relationships'),
               xLabels =colnames(design), yLabels = names(module), ySymbols = names(module),
               colorLabels = FALSE, colors = greenWhiteRed(50), cex.text = 0.7, zlim = c(-1,1),
               textMatrix = textMatrix, setStdMargins = FALSE)
dev.off()
getwd()


##############################################################################################
table(mergedColors)

gene_module <- data.frame(gene_name = colnames(gene), module = mergedColors, stringsAsFactors = FALSE)
head(gene_module)

gene_module_select <- subset(gene_module, module == 'midnightblue')$gene_name

gene_select <- t(gene[,gene_module_select])

tom_select <- tom_sim[gene_module_select,gene_module_select]


#write.table(gene_select, 'gene_select.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.table(tom_select, 'tom_select.txt', sep = '\t', col.names = NA, quote = FALSE)


##################################################################################################

gene_midnightblue <- gene[ ,gene_module_select]

geneModuleMembership <- signedKME(gene_midnightblue, module['MEmidnightblue'], outputColumnName = 'MM')
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(module)))

geneTraitSignificance <- as.data.frame(cor(gene_midnightblue, design[,"Luminal"], use = 'p'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(design)))

geneModuleMembership[geneModuleMembership<0.5 | MMPvalue>0.01] <- 0
geneTraitSignificance[geneTraitSignificance<0.5 | MMPvalue>0.01] <- 0

select <- cbind(geneModuleMembership, geneTraitSignificance)
select <- subset(select, geneModuleMembership>=0.5 & geneTraitSignificance>=0.5)
head(select)

plotNetworkHeatmap(gene,
                   plotGenes =sample(colnames(gene),20)  ,
                   networkType = 'unsigned',
                   useTOM = TRUE,
                   power = powers,
                   main = 'unsigned correlations')




###################################################################################3
#网络可视化
dir.create('cytoscape', recursive = TRUE)

for (mod in 1:nrow(table(mergedColors))) {
  modules <- names(table(mergedColors))[mod]
  probes <- colnames(gene)
  inModule <- (mergedColors == modules)
  modProbes <- probes[inModule]
  modGenes <- modProbes
  modtom_sim <- tom_sim[inModule, inModule]
  dimnames(modtom_sim) <- list(modProbes, modProbes)
  outEdge <- paste('cytoscape/', modules, '.edge_list.txt',sep = '')
  outNode <- paste('cytoscape/', modules, '.node_list.txt', sep = '')
  exportNetworkToCytoscape(modtom_sim,
                           edgeFile = outEdge,
                           nodeFile = outNode,
                           weighted = TRUE,
                           threshold = 0.3,  #该参数可控制输出的边数量，过滤低权重的边
                           nodeNames = modProbes,
                           altNodeNames = modGenes,
                           nodeAttr = mergedColors[inModule])
}






