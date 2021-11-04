###Set up working directory
setwd("//Users/wenying/Desktop/Newest_version_of_manuscript/Resemblance")
#-------------------------------------------------------------------------------#
#Prepare libraryies
library(tidyverse)
library(readxl)
library(writexl)
#-------------------------------------------------------------------------------#
#Read and combine files
PDCPath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/All Primary cell lines vs blood_WES_DNA Link_20191211.xlsx"
PDC <- read_excel(PDCPath, sheet = 1)
TumorPath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/All Fresh tumors vs blood_WES_DNA Link_20191216.xlsx"
Tumor <- read_excel(TumorPath, sheet = 1)

Select = c("ID","Variant_Classification","Hugo_Symbol","Start_position","AF")
PDC <- subset.data.frame(PDC, select = Select)
Tumor <- subset.data.frame(Tumor, select = Select)
All <- rbind.data.frame(PDC,Tumor) #Combine dataframe
All$Hugo <- str_remove(All$Hugo_Symbol,"`") #remove ` in Hugo symbol
#-------------------------------------------------------------------------------#
#Get nonsynonomous data
Type <- unique(as.vector(All$Variant_Classification))
Type
Nonsyn <- Type[c(2,3,5,6,9,10,11,13,14,16,17,20)]
Nonsyn
Allnonsyn <- All %>% filter(Variant_Classification %in% Nonsyn)
#-------------------------------------------------------------------------------#
#QM13
ID = "QM13"
IDC = "QM13-P32"
IDT = "QM13T"
Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list(Same, TumorUnique, LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#QM16
ID = "QM16"
IDC = "QM16-P32"
IDT = "QM16T"

Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#QM17
ID = "QM17"
IDC = "QM17-P31"
IDT = "QM17T"

Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#QM28
ID = "QM28"
IDC = "QM28-P32"
IDT = "QM28T"

Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#QM29
ID = "QM29"
IDC = "QM29-P31"
IDT = "QM29T"

Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#T63
ID = "T63"
IDC = "T63-P32"
IDT = "T63T"

Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#T76
ID = "T76"
IDC = "T76-P20"
IDT = "T76T"

Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#T77
ID = "T77"
IDC = "T77-P26"
IDT = "T77T"


Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)
#-------------------------------------------------------------------------------#
#QM43
ID = "QM43"
IDC = "QM43-P24"
IDT = "QM43T"


Filename = paste(ID,".xlsx")
Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)

Same <- Line %>% filter(Start_position %in% Tu$Start_position) #Tumor and line shared mutations
TumorUnique <- Tu %>% filter(!Start_position %in% Same$Start_position)
LineUnique <- Line %>% filter(!Start_position %in% Same$Start_position)
x <- list("same" = Same, "tumor" = TumorUnique, "PDC" = LineUnique)
write_xlsx(x, path = Filename, col_names = TRUE, format_headers = TRUE)


#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#QM25 early
ID = "QM25_Pd10"
IDAF = "QM25AF_Pd10"
IDC = "QM25-P4"
IDT = "QM25T"

Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(sPDC = Line$Start_position, Tumor= Tu$Start_position)
venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF <- Line %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(sPDC = LineAF$Start_position, Tumor= TuAF$Start_position)
venn.diagram(x, IDAF, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#QM25 middle
ID = "QM25_Pd60"
IDAF = "QM25AF_Pd60"
IDC = "QM25-P32"
IDT = "QM25T"

Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(sPDC = Line$Start_position, Tumor= Tu$Start_position)
venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF <- Line %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(sPDC = LineAF$Start_position, Tumor= TuAF$Start_position)
venn.diagram(x, IDAF, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#QM25 late
ID = "QM25_Pd200"
IDAF = "QM25AF_Pd200"
IDC = "QM25-P82"
IDT = "QM25T"

Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(sPDC = Line$Start_position, Tumor= Tu$Start_position)
venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF <- Line %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(sPDC = LineAF$Start_position, Tumor= TuAF$Start_position)
venn.diagram(x, IDAF, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#T48 early
ID = "T48_Pd10"
IDAF = "T48AF_Pd10"
IDC = "T48-P8"
IDT = "T48T"

Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(sPDC = Line$Start_position, Tumor= Tu$Start_position)
venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF <- Line %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(sPDC = LineAF$Start_position, Tumor= TuAF$Start_position)
venn.diagram(x, IDAF, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#T48 middle
ID = "T48_Pd60"
IDAF = "T48AF_Pd60"
IDC = "T48-P33"
IDT = "T48T"

Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(sPDC = Line$Start_position, Tumor= Tu$Start_position)
venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF <- Line %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(sPDC = LineAF$Start_position, Tumor= TuAF$Start_position)
venn.diagram(x, IDAF, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#T48 late
ID = "T48_Pd200"
IDAF = "T48AF_Pd200"
IDC = "T48-P95"
IDT = "T48T"

Line <- Allnonsyn %>% filter(ID == IDC)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(sPDC = Line$Start_position, Tumor= Tu$Start_position)
venn.diagram(x, ID, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF <- Line %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(sPDC = LineAF$Start_position, Tumor= TuAF$Start_position)
venn.diagram(x, IDAF, fill = c("blue", "red"), cat.pos = c(320,30), height = 900, width = 900, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#T48 three passage
ID = "T48"
Filename = "T483P"
FilenameAF = "T483PAF"
IDC8 = "T48-P8"
IDC33 = "T48-P33"
IDC95 = "T48-P95"

Line8 <- Allnonsyn %>% filter(ID == IDC8)
Line33 <- Allnonsyn %>% filter(ID == IDC33)
Line95 <- Allnonsyn %>% filter(ID == IDC95)
x <- list(Pd10 = Line8$Start_position, Pd60 = Line33$Start_position, Pd200 = Line95$Start_position)
venn.diagram(x, Filename, fill = c("blue","yellow","green"), cat.pos = c(310,50,180), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF8 <- Line8 %>% filter(AF >= 0.1)
LineAF33 <- Line33 %>% filter(AF >= 0.1)
LineAF95 <- Line95 %>% filter(AF >= 0.1)
x <- list(Pd10 = LineAF8$Start_position, Pd60 = LineAF33$Start_position, Pd200 = LineAF95$Start_position)
venn.diagram(x, FilenameAF, fill = c("blue","yellow","green"), cat.pos = c(310,50,180), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#T48 three passage and tumor
ID = "T48"
Filename = "T48Tumor3P"
FilenameAF = "T48Tumor3PAF"
IDC8 = "T48-P8"
IDC33 = "T48-P33"
IDC95 = "T48-P95"
IDT = "T48T"

Line8 <- Allnonsyn %>% filter(ID == IDC8)
Line33 <- Allnonsyn %>% filter(ID == IDC33)
Line95 <- Allnonsyn %>% filter(ID == IDC95)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(Pd10 = Line8$Start_position, Tumor= Tu$Start_position, Pd60 = Line33$Start_position, Pd200 = Line95$Start_position)
venn.diagram(x, Filename, fill = c("blue", "red","yellow","green"), cat.pos = c(0,0,0,0), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF8 <- Line8 %>% filter(AF >= 0.1)
LineAF33 <- Line33 %>% filter(AF >= 0.1)
LineAF95 <- Line95 %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(Pd10 = LineAF8$Start_position, Tumor= TuAF$Start_position, Pd60 = LineAF33$Start_position, Pd200 = LineAF95$Start_position)
venn.diagram(x, FilenameAF, fill = c("blue", "red","yellow","green"), cat.pos = c(0,0,0,0), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#QM25 three passage
ID = "QM25"
Filename = "QM253P"
FilenameAF = "QM253PAF"
IDC8 = "QM25-P4"
IDC33 = "QM25-P32"
IDC95 = "QM25-P82"

Line8 <- Allnonsyn %>% filter(ID == IDC8)
Line33 <- Allnonsyn %>% filter(ID == IDC33)
Line95 <- Allnonsyn %>% filter(ID == IDC95)
x <- list(Pd10 = Line8$Start_position, Pd60 = Line33$Start_position, Pd200 = Line95$Start_position)
venn.diagram(x, Filename, fill = c("blue","yellow","green"), cat.pos = c(310,50,180), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF8 <- Line8 %>% filter(AF >= 0.1)
LineAF33 <- Line33 %>% filter(AF >= 0.1)
LineAF95 <- Line95 %>% filter(AF >= 0.1)
x <- list(Pd10 = LineAF8$Start_position, Pd60 = LineAF33$Start_position, Pd200 = LineAF95$Start_position)
venn.diagram(x, FilenameAF, fill = c("blue","yellow","green"), cat.pos = c(310,50,180), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)
#-------------------------------------------------------------------------------#
#QM25 three passage and tumor
ID = "QM25"
Filename = "QM253PTumor"
FilenameAF = "QM253PAFTumor"
IDC8 = "QM25-P4"
IDC33 = "QM25-P32"
IDC95 = "QM25-P82"
IDT = "QM25T"

Line8 <- Allnonsyn %>% filter(ID == IDC8)
Line33 <- Allnonsyn %>% filter(ID == IDC33)
Line95 <- Allnonsyn %>% filter(ID == IDC95)
Tu <- Allnonsyn %>% filter(ID == IDT)
x <- list(Pd10 = Line8$Start_position, Tumor= Tu$Start_position, Pd60 = Line33$Start_position, Pd200 = Line95$Start_position)
venn.diagram(x, Filename, fill = c("blue", "red","yellow","green"), cat.pos = c(0,0,0,0), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

LineAF8 <- Line8 %>% filter(AF >= 0.1)
LineAF33 <- Line33 %>% filter(AF >= 0.1)
LineAF95 <- Line95 %>% filter(AF >= 0.1)
TuAF <- Tu %>% filter(AF >= 0.1)
x <- list(Pd10 = LineAF8$Start_position, Tumor= TuAF$Start_position, Pd60 = LineAF33$Start_position, Pd200 = LineAF95$Start_position)
venn.diagram(x, FilenameAF, fill = c("blue", "red","yellow","green"), cat.pos = c(0,0,0,0), height = 1000, width = 1000, resolution = 200, imagetype = "tiff", units = "px", compression ="lzw", na = "stop", main = ID, sub = NULL, main.pos= c(0.5, 1), main.fontface = "plain",main.fontfamily = "serif", main.col = "black",main.cex = 1.5, main.just = c(0.5, 1), category.names = names(x), force.unique =TRUE, print.mode = "raw", sigdigs = 3, direct.area =FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, lower.tail = TRUE)

