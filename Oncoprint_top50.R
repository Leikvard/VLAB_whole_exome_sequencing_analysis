###Set up working directory
setwd("//Users/wenying/Desktop/Newest_version_of_manuscript/Oncoprint/Common")
#-------------------------------------------------------------------------------#
#Prepare libraryies
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
#-------------------------------------------------------------------------------#
#column_title = "Significantly mutated genes in HNSCC"
#Filename1 <- "PI3K .csv"
#-------------------------------------------------------------------------------#
#Read files
GenePath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/Oncoprint/Common/Top50.xlsx"
GeneList <- read_excel(GenePath)
GeneList <- as.vector(GeneList$V1)
GeneList

PDCPath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/All Primary cell lines vs blood_WES_DNA Link_20191211.xlsx"
PDC <- read_excel(PDCPath, sheet = 1)
TumorPath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/All Fresh tumors vs blood_WES_DNA Link_20191216.xlsx"
Tumor <- read_excel(TumorPath, sheet = 1)

All <- rbind.data.frame(PDC,Tumor) #Combine dataframe
All$Hugo <- str_remove(All$Hugo_Symbol,"`") #remove ` in Hugo symbol
#-------------------------------------------------------------------------------#
#Prepare nonsynonymous list
TypeDefPath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/NONSYN SYN RNA mutations DNA link consensus_VL_19Dec2019.xlsx"
TypeDef <- read_excel(TypeDefPath, sheet = 2)
NonsynType <- TypeDef$`Non-syn`
Nonsyn <- All %>% filter(Variant_Classification %in% NonsynType)
NSDat <- Nonsyn %>% filter(Hugo %in% GeneList) #Filter with targetable genes in dataset
#-------------------------------------------------------------------------------#
#Prepare Tumor ID list
IDPath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/Label_name.xlsx"
ID <- read_excel(IDPath, sheet = 1)
IDList <- as.vector(ID$V1)
IDLabel <- as.vector(ID$V2)
#IDLabel[is.na(IDLabel)] <- ""
#-------------------------------------------------------------------------------#
#Prepare output table
Output <- matrix(, nrow = length(GeneList), ncol = length(IDList))
Output <- as.data.frame(Output)
rownames(Output) <- GeneList
colnames(Output) <- IDList
Output #check Output before processing
#-------------------------------------------------------------------------------#
#Read MUT into output
MutID <- unique(NSDat$ID)
for (i in c(1:length(MutID))){
  Sample <- MutID[i]
  df <- NSDat %>% filter(ID == Sample)
  for (j in c(1 : nrow(df))){
    Gene <- as.vector(df$Hugo)[j]
    Type <- as.vector(df$Variant_Classification)[j]
    Output[Gene,Sample] = paste(Output[Gene,Sample],";",Type)
  }
}
Filename <- "Output20200113_2.csv"
write.csv(Output, file = Filename, sep = "\t", na = "")

#-------------------------------------------------------------------------------#
#####ONCOPRINT######
#-------------------------------------------------------------------------------#
#Oncoprint preparations
df <- as.matrix(Output20200113_2)
df[is.na(df)] = ""
rownames(df) = df[, 1]
df = df[, -1]

#Fill empty column
EmptyCol <- letters[1:12]
for (i in c(1:length(EmptyCol))){
  x <- EmptyCol[i]
  for (j in c(1:nrow(Output))){
    df[j,x] <- "white"
  }
}

#The type name in col, Type and alter_fun variable has to be the same as in the file
col = c("Missense" = "green4", "Nonsense" = "black","Inframe Indel" = "brown","Frameshift Indel" = "chocolate", "Others" = "coral", "white" = "white")
Type = c("Missense","Nonsense","Inframe Indel","Frameshift Indel","Others","white")
TypeLabel = c('missense','truncating',"inframe indel","frameshift indel", "others","")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # lightgreen
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.9, 
              gp = gpar(fill = col["Missense"], col = NA))
  },
  # black
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = col["Nonsense"], col = NA))
  },
  # brown
  `Inframe Indel` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = col["Inframe Indel"], col = NA))
  },
  # chocolate
  `Frameshift Indel` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = col["Frameshift Indel"], col = NA))
  },
  # chocolate
  Others = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = col["Others"], col = NA))
  },
  # chocolate
  white = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*1, 
              gp = gpar(fill = col["white"], col = NA))
  }
)

heatmap_legend_param = list(title = "Alterations", at = Type, 
                            labels = TypeLabel)

ha = HeatmapAnnotation(
  text = anno_text(IDLabel, location = 1, rot = 90, just = 'right')
)

oncoPrint(df,
          alter_fun = alter_fun, col = col, remove_empty_rows = TRUE,
          bottom_annotation = ha, column_order = IDList, row_names_gp = gpar(fontface = "italic",size = 0.1),
          heatmap_legend_param = heatmap_legend_param)



###Column annotation
anno_df <- read_excel(IDPath, sheet = 2)
df[is.na(df)] <- ""
gender = anno_df[, "Gender"]
yearstobirth = anno_df[, "Age"]
pathologicstage = anno_df[, "Stage"]
ha_ca = HeatmapAnnotation(gender = gender, stage = pathologicstage,
                       #age = anno_points(yearstobirth, ylim = c(0, max(yearstobirth, na.rm = TRUE)), axis = TRUE),
                       col = list(gender = c("male" = "red", "female" = "blue", "NA" = "white"),
                                  stage = c("stage ii" = "#60FF60", "stage iii" = "#B0FFB0", "NA" = "white",
                                            "stage iva" = "#6060FF", "stage ivb" = "#B0B0FF",
                                            "stage iv" = "#FFFF00")),
                       annotation_height = unit(c(5, 15), "mm"),
                       annotation_legend_param = list(gender = list(title = "Gender"),
                                                      stage = list(title = "Stage"))
)
pathologicstage
gender




