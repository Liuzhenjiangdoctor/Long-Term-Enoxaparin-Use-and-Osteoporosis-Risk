rm(list=ls())
library(tidyverse)
setwd('C:/Users/LZJ/Desktop/File')
load('GSE35958.rda')
dif <- read.table('DEG2.txt',header = TRUE, sep = "\t",row.names = 1)
difs <- rownames(dif)[dif$change!='NOT']
dat <- t(exp[difs,])
rt <- rt
rt$group <- factor(rt$group, levels = c("control", "pSS"))
identical(rownames(rt),rownames(dat))
library(glmnet
set.seed(123)
x <- as.matrix(dat)
y <- rt$group
fit <- glmnet(x,y,family='binomial', alpha=1)
pdf(file="lambda.pdf", height= 4.5, width= 4.5)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit = cv.glmnet(x, y, alpha = 1, family="binomial",nfolds = 10)
pdf(file="deviance.pdf", height= 4.5, width= 4.5)
plot(cvfit)
dev.off()
coef <- coef(fit,s=cvfit$lambda.min)
index <- which(coef!=0)
actCoef <- coef[index]
lassogenes <- row.names(coef)[index]
lassogenes
geneCoef <- cbind(Gene=lassogenes,Coef=actCoef)
geneCoef <- as.data.frame(geneCoef[-1,])
order <- arrange(geneCoef,Coef)
order$Coef <- as.numeric(order$Coef)
par(mar = c(5, 8, 4, 2)) 
barplot(order$Coef,names.arg = order$Gene,
        horiz=TRUE,las=1,
        col=ifelse(order$Coef > 0,"tomato","steelblue"),
        xlab="Coefficient Value",
        main="Biomarker Impact Direction",
        cex.names = 0.9,cex.lab = 1,cex.main = 1)