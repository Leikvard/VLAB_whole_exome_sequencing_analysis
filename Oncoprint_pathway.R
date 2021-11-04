###Set up working directory
setwd("//Users/wenying/Desktop/Newest_version_of_manuscript/Oncoprint")
#-------------------------------------------------------------------------------#
#Prepare libraryies
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
#-------------------------------------------------------------------------------#
#Read files
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
#-------------------------------------------------------------------------------#
#Prepare Tumor ID list
IDPath <- "//Users/wenying/Desktop/Newest_version_of_manuscript/Label_name.xlsx"
ID <- read_excel(IDPath, sheet = 1)
IDList <- as.vector(ID$V1)
IDLabel <- as.vector(ID$V2)
IDLabel[is.na(IDLabel)] <- ""
#-------------------------------------------------------------------------------#
#Prepare gene list and data for gene list
GeneList <- read_excel("//Users/wenying/Desktop/Newest_version_of_manuscript/gene_list/all_vlab.xlsx", sheet = 4)
GeneList <- as.vector(GeneList$NFkB)
GeneList
NSDat <- Nonsyn %>% filter(Hugo %in% GeneList) #Filter with targetable genes in dataset
#-------------------------------------------------------------------------------#
#Prepare output table
Output <- matrix(, nrow = length(GeneList), ncol = length(IDList))
Output <- as.data.frame(Output)
rownames(Output) <- GeneList
colnames(Output) <- IDList
#-------------------------------------------------------------------------------#
#Read MUT into output
MutID <- unique(NSDat$ID)
for (i in c(1:length(MutID))){
  Sample <- MutID[i]
  df <- NSDat %>% filter(ID == Sample)
  for (j in c(1 : nrow(df))){
    Gene <- as.vector(df$Hugo)[j]
    Type <- as.vector(df$Variant_Classification)[j]
    if (Type == "Missense_Mutation"){
      Output[Gene,Sample] = paste(Output[Gene,Sample],";","missense")
    }
    else if (Type == "Nonsense_Mutation"){
      Output[Gene,Sample] = paste(Output[Gene,Sample],";","truncating")
    }
    else if (Type == "Frame_Shift_Ins" | Type == "Frame_Shift_Del" | Type == "In_Frame_Ins" | Type == "In_Frame_Del"){
      Output[Gene,Sample] = paste(Output[Gene,Sample],";","indel")
    }
    else {
      Output[Gene,Sample] = paste(Output[Gene,Sample],";","others")
    }
  }
}
#-------------------------------------------------------------------------------#
#####ONCOPRINT######
#-------------------------------------------------------------------------------#
#Oncoprint preparations
df <- as.matrix(Output)
df <- df[rowSums(is.na(df)) != ncol(df), ]
df[is.na(df)] = ""
#Fill empty column
EmptyCol <- letters[1:12]
for (i in c(1:length(EmptyCol))){
  x <- EmptyCol[i]
  for (j in c(1:nrow(df))){
    df[j,x] <- "white"
  }
}

#The type name in col, Type and alter_fun variable has to be the same as in the file
col = c("missense" = "green4", "truncating" = "black","indel" = "coral2", "others" = "cornflowerblue", "white" = "white")
Type = c("missense","truncating","indel","others","white")
TypeLabel = c("missense","truncating","indel", "others","")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # lightgreen
  missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.9, 
              gp = gpar(fill = col["missense"], col = NA))
  },
  # black
  truncating = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = col["nonsense"], col = NA))
  },
  # brown
  indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = col["indel"], col = NA))
  },
  # chocolate
  others = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = col["others"], col = NA))
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




