Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(GEOquery)
library(dplyr)
library(tidyverse)

gset = getGEO(GEO='GSE12417', destdir=".",getGPL = F)
##save(gset,file = "gse12417.rdata")
##load("gse12417.rdata")
e2=gset[[2]]
exp2=e2@assayData[["exprs"]]

phe=pData(e2)
pdata1=e2@phenoData@data
#########################################################################################################
library(data.table)
anno=fread("GPL96-57554.txt",sep = "\t",header = T,data.table = F)
colnames(anno)
gene=anno[,c("ID","Gene Symbol")]
exp1=as.data.frame(exp2)
exp.anno=merge(x=gene,y=exp1,by.x=1,by.y=0)
x=gene$`Gene Symbol`

a1=strsplit(x,split = " /// ",fixed = T)
gene.all = sapply(a1,function(x){x[1]})
exp.anno$`Gene Symbol`=gene.all
###########################################################################################

exp1=exp.anno
colnames(exp1)[2]="gene.all"

exp2=distinct(exp1,gene.all,.keep_all = T)

exp3=na.omit(exp2)
rownames(exp3)=exp3$gene.all
View(exp3)
exp4=exp3[,-c(1,2)]

######################
####生存曲线
phe=e2@phenoData@data
s=phe$characteristics_ch1
s[1:5]
s1=strsplit(s,split = "; ",fixed = T)
s1[[1]][1:5]
os.time = sapply(s1,function(x){x[3]})
status= sapply(s1,function(x){x[4]})
os.time[1:5]
os.time1=strsplit(os.time,split = " ",fixed = T)
os.time2= sapply(os.time1,function(x){x[3]})
os.time2=as.numeric(os.time2)
status[1:5]
status1=strsplit(status,split = ": ",fixed = T)
status2= sapply(status1,function(x){x[2]})
status2=as.numeric(status2)
s2=data.frame(os.time2,status2)
rownames(s2)=rownames(phe)
######################################################

exp4.e=exp4["EZH2",]

exp4.e <- data.frame(t(exp4.e))

###s2.exp <- merge(s2,exp4.e, by.x=0, by.y=0)
s2.exp <- merge(s2,exp4.e, by=0)


##############################################################################
###四分位生存曲线
quantile(s2.exp$EZH2)
s2.exp$zz=ifelse(s2.exp$EZH2>9.53,'high', 
                 ifelse( s2.exp$EZH2>9.23 & s2.exp$EZH2<9.53,'h.stable', 
                         ifelse( s2.exp$EZH2>8.88& s2.exp$EZH2<9.23,'l.stable','down') )
)
############################################################################
##中位数分组
median(s2.exp$EZH2)
s2.exp$group = ifelse(     s2.exp$EZH2 > 9.235157 ,'high','low')


s2.exp$EZH2 = ifelse(     s2.exp$EZH2 > 9.235157           ,'H','low')

################## 

library(survival)
library(survminer)

############################################
Surv(os.time2,status2)

x1=survfit(Surv(os.time2, status2)~group, data=s2.exp)
summary(x1)

plot(x1)

plot(x1,xlab="Time(Days)",ylab="Survival probability",
     col=c("blue","red"),lty=2:3,lwd=2) 
legend("topright",c("high","low"),
       col=c("blue","red"),lty=2:3,lwd=2,cex=1)



################################################################################

ggsurvplot(x1, data = s2.exp, pval = T,
           ylab="Fraction Without Recurrence",xlab="Time (days)",
           linetype=c("dashed","solid"),palette="jco",risk.table=T,
           pval.coord = c(90, 0.9))



ggsurvplot(x1, data = s2.exp, pval = T,ylab="Survival probability",conf.int=T,
           xlab="Time(Days)",linetype=c("dashed","solid"),palette="lancet",risk.table=T,pval.coord = c(90, 0.9))


ggsurvplot(fit)

ggsurvplot(survfit(Surv(os.time2, status2)~zz, data=s2.exp), conf.int=T, pval=TRUE)

ggsurvplot(survfit(Surv(os.time2, status2)~zz, data=s2.exp), conf.int=F, pval=T)


ggsurvplot(x1,   
           main = "Survival curve", # 添加标题
           font.main = c(16, "bold", "darkblue"), # 设置标题字体大小、格式和颜色
           font.x = c(14, "bold.italic", "red"), # 设置x轴字体大小、格式和颜色
           font.y = c(14, "bold.italic", "darkred"), # 设置y轴字体大小、格式和颜色
           font.tickslab = c(12, "plain", "darkgreen")) # 设置坐标轴刻度字体大小、格式和颜色
######################################################

ggsurvplot(x1, 
           surv.median.line = "hv", # 添加中位数生存时间线
           
           # Change legends: title & labels
           legend.title = "EZH2", # 设置图例标题
           legend.labs = c("high", "low"), # 指定图例分组标签
           
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = TRUE, # 设置添加置信区间
           
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_gray(), # 设置风险表的主题
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#FFFF00", "#0000AA"), # 设置颜色画板
           ggtheme = theme_classic() # Change ggplot2 theme
)
?ggsurvplot

################################
####十六进制颜色

fit <- survfit(Surv(os.time2, status2)~EZH2, data=s2.exp)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#FF0000", "#5599FF")
)

################################
####十六进制颜色

fit <- survfit(Surv(os.time2, status2)~zz, data=s2.exp)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw()  # Change ggplot2 theme
           
)

####################################################

# 使用ggsurvplot_facet()函数绘制分面生存曲线
s2.exp$sex=c(rep("m",100),rep("f",63))
fit=survfit(Surv(os.time2, status2)~zz, data=s2.exp)
ggsurvplot_facet(fit, s2.exp, 
                 facet.by = "sex", # 设置分面变量
                 palette = "jco", # 设置颜色画板
                 pval = TRUE) # 添加pvalue值
#########################换个配色
ggsurvplot_facet(fit, s2.exp,  
                 facet.by = c("sex"),
                 palette = "npg", 
                 pval = TRUE,
                 surv.median.line = "hv",  # 增加中位生存时间
                 conf.int = TRUE) # 增加置信区间)

###################################
# 拟合多个分组变量
fit2 <- survfit( Surv(os.time2, status2) ~ sex + zz, data = s2.exp )
fit2

ggsurvplot(fit2)

#####绘制累计风险曲线
ggsurvplot(fit, data = s2.exp, 
           conf.int = TRUE, # 增加置信区间
           fun = "cumhaz") # 绘制累计风险曲线

#############添加总患者生存曲线

ggsurvplot(fit, # 创建的拟合对象
           data = s2.exp,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           add.all = TRUE) # 添加总患者生存曲线




ggsurvplot(fit2,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           
           fun = "event")
res.sum <- surv_summary(fit)
##想看fit每个时间点的结果用 surv_summary（）函数
ggsurvplot(res.sum,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")
)

