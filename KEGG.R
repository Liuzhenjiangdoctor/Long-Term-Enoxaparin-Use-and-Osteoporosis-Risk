#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")
#install.packages("ggpubr")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("ComplexHeatmap")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
pvalueFilter=0.05        
p.adjustFilter=0.05      
colorSel="p.adjust"
if(p.adjustFilter>0.05){
  colorSel="pvalue"
}
setwd("C:\\Users\\lexb\\Desktop\\Toxicity\\17.KEGG")      
rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)    
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<p.adjustFilter),]
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)
showNum=30     
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}
pdf(file="barplot.pdf", width=8, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()
pdf(file="bubble.pdf", width=8, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()
pdf(file="Lollipop.pdf", width=9.5, height=6)
ggdotchart(KEGG, x="Description", y="Count", color = "category",group = "category", xlab="",
           palette = "aaas",     
           legend = "right",     
           sorting = "descending",   
           add = "segments",    
           rotate = TRUE,       ʾ
           dot.size = 7,       
           label = round(KEGG$Count),   
           font.label = list(color="white",size=12, vjust=0.5),   
           ggtheme = theme_pubr())
dev.off()

