##清理R语言的所有对象##
rm(list = ls(all.names = TRUE))
##遇到字符串之后，不将其转换为factors，仍然保留为字符串格式##
options(stringsAsFactors = F)
##矩阵稀疏的时候使用Matrix节约空间##
library(Matrix)
##设置工作路径##
setwd('D:\\SCIwork\\F38KRT\\s2')

##调整表达数据格式，增加样品名作为新的一列##
data <- read.csv('cdata.csv', header = T, row.names = 1)

data <- as.data.frame(t(data))

data[1:4,1:4]

data$sample <- rownames(data)

data$sample <- chartr(old='.', new='-', x=data$sample)

##调整样品信息，合并表达数据和样品信息##
setwd('D:\\SCIwork\\F38KRT\\s3')

group <- read.csv('group2.csv', header = T)

names(group)[1] <- 'sample'

group$sample <- chartr(old='.', new='-', x=group$sample)

group <- subset(group, select=c("sample", "group"))

group$subtype <- group$group

group$group  <- NULL

dt <- merge(group, data, by='sample')

dt[1:4,1:4]

dt$sample <- NULL

table(dt$subtype)
##保存一个dt_total,方便后面使用##
dt_total <- dt

#===========================================================================
##进行T检验之前的分组调整，一组为实验组，其余所有组变为control##
#===========================================================================
dt$subtype <- NULL
dt <- as.data.frame(t(dt))
dt <- rowSums(dt> 0) >= floor(0.001*ncol(dt))
dt <- as.data.frame(t(dt))
dt$subtype <- dt_total$subtype

dt$subtype <- ifelse(dt$subtype == 'Subtype1', 'Exp', 'Con')

table(dt$subtype)

dt <- dt[order(dt$subtype), ]

dt[1:4,1:4]

dt_Con <- subset(dt, dt$subtype == 'Con')

dt_Con[1:4,1:4]

dt_Exp <- subset(dt, dt$subtype == 'Exp')

dt_Exp[1:4,1:4]



##subtype信息变为subtype+样本名##
dt_Con$subtype <- paste0(dt_Con$subtype, rownames(dt_Con))

rownames(dt_Con) <- dt_Con$subtype

dt_Con$subtype <- NULL
##数据变为每个基因为行，样品为列##
dt_Con <- as.data.frame(t(dt_Con))


dt_Exp$subtype <- paste0(dt_Exp$subtype, rownames(dt_Exp))

rownames(dt_Exp) <- dt_Exp$subtype

dt_Exp$subtype <- NULL

dt_Exp <- as.data.frame(t(dt_Exp))


##设置Pvalue向量##
Pvalue<-c(rep(0,nrow(dt_Con)))
##设置差异倍数向量##
log2_FC<-c(rep(0,nrow(dt_Con)))
##用for循环进行每一行基因的t检验，得到pvalue和差异倍数##
for(i in 1:nrow(dt_Con)){
  
  y=t.test(as.numeric(dt_Con[i,]),as.numeric(dt_Exp[i,]))
  Pvalue[i] <- y$p.value
  log2_FC[i] <-log2(mean(as.numeric(dt_Exp[i,]))/(mean(as.numeric(dt_Con[i,]))))
  
}


library(dplyr)

library(tidyr)

library(tibble)

# 对pvalue进行FDR校正
fdr=p.adjust(Pvalue, "BH") 
# 在原文件后面加入log2FC，p value和FDR,共3列；
out<- as.data.frame(cbind(log2_FC,Pvalue,fdr))
out$gene <- rownames(dt_Con)
# out <- out %>%
out <- filter(out,log2_FC > 0)
out <-out [order(out$Pvalue),] 
out <- out[1:floor(0.1*nrow(out)),]

setwd('D:\\SCIwork\\F38KRT\\s5')

write.csv(out, file = 'out_S1.csv')
