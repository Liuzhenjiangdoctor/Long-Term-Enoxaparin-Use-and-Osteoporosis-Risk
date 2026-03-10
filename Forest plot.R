setwd("C:/Users/LZJ/Desktop/File")
getwd()
library(grid)
library(forestploter)
library(openxlsx)
dt <- read.csv("xlxs.csv",sep=',',header = T)
dt<-dt[,1:6]
dt$Variables <- ifelse(is.na(dt$or), 
                       dt$Variables,
                       paste0("   ", dt$Variables))
dt$` ` <- paste(rep(" ", 20), collapse = " ")
dt$'OR(95%CI)'<-ifelse(is.na(dt$or),"",
                       sprintf('%.2f(%.2f to %.2f)',
                               dt$or,dt$or_lci95,dt$or_uci95))
library(grid) # Required for gpar() function
tm <- forest_theme(base_size = 10,
                   ### CI settings (unchanged)
                   ci_pch = 20,
                   ci_col = "#4575b4",
                   ci_lty = 1,
                   ci_lwd = 2.3,
                   ci_Theight = 0.2,
                   ### Reference line (Updated to refline_gp)
                   # Old: refline_lwd, refline_lty, refline_col
                   refline_gp = gpar(lwd = 1.5, lty = "dashed", col = "red"),
                   
                   ### Summary diamonds (unchanged)
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   
                   ### Footnote (Updated to footnote_gp)
                   # Old: footnote_cex, footnote_fontface, footnote_col
                   footnote_gp = gpar(cex = 1.1, fontface = "italic", col = "blue")
)
p <- forest(dt[,c(1,7:8, 2:3)],
            est = dt$or,
            lower = dt$or_lci95,
            upper = dt$or_uci95,
            sizes = 0.6,
            ci_column = 2,
            ref_line = 1,
            xlim = c(0,4),
            ticks_at = c(0,1,2,3,4),
            arrow_lab = c('protective factor','risk factor'),
            theme = tm)
print(p)