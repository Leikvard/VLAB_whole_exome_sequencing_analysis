###Set up initial values
setwd("//Users/wenying/Desktop")
#-------------------------------------------------------------------------------#
#Prepare libraryies
library(tidyverse)
library(readxl)
library(VennDiagram)
#-------------------------------------------------------------------------------#
#Read files and construct dataframes
PDCPath <- "//Users/wenying/Desktop/All Primary cell lines vs blood_WES_DNA Link_20191211.xlsx"
PDC <- read_excel(PDCPath, sheet = 1)
TumorPath <- "//Users/wenying/Desktop/All Fresh tumors vs blood_WES_DNA Link_20191216.xlsx"
Tumor <- read_excel(TumorPath, sheet = 1)
#-------------------------------------------------------------------------------#
#Draw 2 density plots and Venn plot, change ID, IDC and IDT
ID = "T77"
IDC = "T77-P26"
IDT = "T77T"
LineData <- PDC %>% filter(Cell_line == IDC)
TuData <- Tumor %>% filter(Tumor == IDT)
LineAF <- LineData$AF
TuAF <- TuData$AF
#x <- list(sPDC = LineData$Start_position, Patient_tumor = TuData$Start_position)
#venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 3000, width = 3000, resolution = 500, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
DL <- density(LineAF)
DT <- density(TuAF)
plot(DL, col = "blue", main = ID, xlim = c(0, max(DL$x,DT$x)), ylim = c(0, max(DL$y,DT$y)))
lines(DT, col = "red")
#legend(0.5 * max(DL$x,DT$x), 0.5 * max(DL$y,DT$y), legend=c("sPDC", "Patient tumor"),col=c("blue", "red"), lty=1:1, cex=0.8)
legend(0.7, 2,legend=c("sPDC", "Patient tumor"),col=c("blue", "red"), lty=1, cex=0.8)
#-------------------------------------------------------------------------------#
#Draw 4 density plots and Venn plot, change ID, IDC and IDT
ID = "T48"
IDC1 = "T48-P8"
IDC2 = "T48-P33"
IDC3 = "T48-P95"
IDT = "T48T"
LineData1 <- PDC %>% filter(Cell_line == IDC1)
LineData2 <- PDC %>% filter(Cell_line == IDC2)
LineData3 <- PDC %>% filter(Cell_line == IDC3)
TuData <- Tumor %>% filter(Tumor == IDT)
LineAF1 <- LineData1$AF
LineAF2 <- LineData2$AF
LineAF3 <- LineData3$AF
TuAF <- TuData$AF
#x <- list(sPDC = LineData$Start_position, Patient_tumor = TuData$Start_position)
#venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 3000, width = 3000, resolution = 500, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
DL1 <- density(LineAF1)
DL2 <- density(LineAF2)
DL3 <- density(LineAF3)
DT <- density(TuAF)
plot(DL1, col = "blue", main = ID, xlim = c(0, max(DL$x,DT$x)), ylim = c(0, max(DL$y,DT$y)))
lines(DT, col = "red")
lines(DL2, col = "chocolate4")
lines(DL3, col = "darkgreen")
legend(0.6, 2,legend=c("Early passage", "Middle passage", "Late passage", "Patient tumor"),col=c("blue", "chocolate4","darkgreen", "red"), lty=1, cex=0.8)
