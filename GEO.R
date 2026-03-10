rm(list=ls())
# BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
setwd('C:/Users/LZJ/Desktop/File') 
gset <- getGEO("GSE35958", destdir = ".", AnnotGPL = F, getGPL = F)
class(gset) 
gset[[1]] 
exp <- exprs(gset[[1]])
boxplot(exp, outline=F,las=2)
library(limma)
exp <- normalizeBetweenArrays(exp)
boxplot(exp, outline=F,las=2)
range(exp)
exp <- log2(exp+1)
range(exp)
gpl<-read.table('GPL570-55999.txt',header=TRUE,fill=T,sep='\t',comment.char='#',stringsAsFactors=FALSE,
                quote='')  
View(gpl)
colnames(gpl) 
ids <- gpl[,c('ID','Gene.Symbol')]  
colnames(ids) <- c('probe_id','symbol') 
library(stringr)
ids$symbol <- trimws(str_split(ids$symbol,'///',simplify = T)[,1]) 
ids <- ids[ids$symbol!='',]
ids <- ids[ids$symbol!='---',]
ids <- ids[ids$symbol!='--- ',]
ids <- ids[ids$probe_id %in% rownames(exp),]
ids$probe_id <- as.character(ids$probe_id) 
exp2 <- exp[ids$probe_id,]
table(rownames(exp2)==ids$probe_id)
exp3 <- cbind(ids,exp2)
exp3 <- aggregate( . ~ symbol,exp3,max) 
rownames(exp3) <- exp3$symbol
exp3 <- exp3[,-c(1,2)]
rt1 <- pData(gset[[1]])
write.csv(rt1, file='GSE35958_cli.csv',row.names = T)
colnames(rt1)
group_list <- ifelse(rt1$`source_name_ch1`== 'Control' , "Control","OP")
# group_list <- ifelse(str_detect(rt1$`source_name_ch1`, "Nor"), "Control","Tumor") 
# group_list <- ifelse(str_detect(rt1$`title`, "normal_adjace"), "Normal",ifelse(str_detect(rt1$`title`, "positiveLN_RAIrefractiv"), "P_Tumor","N_Tumor"))
group_list
rt1$group <- group_list
rt2 <- rt1 %>% dplyr::select("group")
rt3 <- rt2 %>% arrange(rt2) 
exp3 <- exp3[,rownames(rt3)]
identical(rownames(rt3),colnames(exp3))
save(rt3,exp3,file = "GSE35958.rda")
rm(list=ls())
load( "GSE35958.rda")
