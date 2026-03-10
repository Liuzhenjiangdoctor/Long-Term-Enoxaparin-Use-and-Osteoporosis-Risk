#install.packages("ggvenn")
library(ggvenn)
setwd("C:\\Users\\lexb\\Desktop\\Toxicity\\09.compound")  
files=list.files(pattern="*.txt$") 
geneList=list()
for(inputFile in files){
  if(inputFile=="compound.txt"){next}
  rt=read.table(inputFile, header=F, sep="\t", check.names=F)      
  geneNames=unlist(strsplit(as.vector(rt[,1]), " "))  
  geneNames=gsub("^ | $","",geneNames)   
  uniqGene=unique(geneNames)           
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#E41A1C","#1E90FF","#FF8C00"),
       set_name_color =c("#E41A1C","#1E90FF","#FF8C00"),
       set_name_size=6, text_size=4.5)
dev.off()
unionGenes=Reduce(union, geneList)
write.table(file="Compound.txt", unionGenes, sep="\t", quote=F, col.names=F, row.names=F)
