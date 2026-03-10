rm(list = ls())
setwd('C:/Users/LZJ/Desktop/File') 
library(limma)
library(pheatmap)
library(tidyverse)
library(ggpubr)
load('GSE35958.rda')
rt <- rt3
exp <- exp3
str(exp) 
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
identical(rownames(rt),colnames(exp))
group_list <- factor(rt$group,levels = c("hMSC-old_elderly","hMSC-osteopo_elderly"))
group_list
design <- model.matrix(~group_list)
design
fit <- lmFit(exp,design)
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG = na.omit(deg) #differently expressed genes
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
logFC_cutoff <- 1
type1 = (DEG$adj.P.Val < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$adj.P.Val < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(DEG$change)
library(ggplot2)
library(cowplot)
library(ggrepel)
DEG <- mutate(DEG,Gene_symbol=row.names(DEG))
UPDEG <- subset(DEG,change=='UP')
UPDEG_5 <- top_n(x = UPDEG,n = -5,wt = P.Value)
DOWNDEG <- subset(DEG,change=='DOWN')
DOWNDEG_5 <- top_n(x = DOWNDEG,n = -5,wt = P.Value)
p <- ggplot(data = DEG,
            aes(x = logFC,
                y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.5, size = 4.5,
             aes(color = change)) +
  ylab("-log10(adj.P.Val)") +
  scale_color_manual(values = c("#1F77B4", "grey", "#FF7F0E")) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5) +
  theme_half_open() +
  geom_text_repel(data = UPDEG_5, aes(label = Gene_symbol), vjust = 1.5, size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50',
                  nudge_y = 0.1,  
                  direction = "y") + 
  geom_text_repel(data = DOWNDEG_5, aes(label = Gene_symbol), vjust = 1.5, size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50',
                  nudge_y = 0.1,  
                  direction = "y") 
p
pdf(file = 'volcano.pdf',width = 10,height = 6)
p
dev.off()
library(pheatmap)
identical(rownames(rt),colnames(exp))
DEG <- read.table("DEG2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
table(DEG$change)
diff <- rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[diff,]
range(exp_diff)
pdf(file = 'heatmap.pdf',width=6,height=20)
pheatmap(exp_diff,
         annotation_col=rt,
         color = colorRampPalette((c("blue","white","red")))(100),
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = F,
         show_rownames = T ,
         border_color = NA,
         scale='row') 
dev.off()
library(tidyverse)
diff =  DEG[DEG$change !="NOT",]
up <- diff %>% top_n(25,logFC)
dw <- diff %>% top_n(-25,logFC)
all <-c(rownames(up),rownames(dw))
exp_diff2 <- exp[all,]
pdf(file = 'heatmap2.pdf',width=6,height=8)
pheatmap(exp_diff2,
         annotation_col=rt,
         color = colorRampPalette((c("blue","white","red")))(100),
         cluster_cols = F,
         cluster_rows = F,
         show_colnames = F,
         show_rownames = T,
         border_color = NA,
         scale='row') 
dev.off()
