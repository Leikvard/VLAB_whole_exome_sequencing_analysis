library(tidyverse)
library(readxl)
CCLEAs <- read.csv(file = "COSMIC_SIGNATURE_CCLE_HNSC_lines_Asia.csv")
rownames(CCLEAs) <- CCLEAs$X
CCLEAs <- CCLEAs[,-1]
Sum <- matrix(0, nrow = nrow(CCLEAs))
for (i in c(1 : nrow(CCLEAs))){
  a = 0
  for (j in c(1 : ncol(CCLEAs))){
    if(CCLEAs[i,j] != 0){
      a = a + 1
    }
  }
  Sum[i,1] = a
}
Sum <- as.data.frame(Sum)
colnames(Sum) <- "CCLEAs"

CCLECa <- read.csv(file = "COSMIC_SIGNATURE_CCLE_HNSC_lines_Caucasia.csv")
rownames(CCLECa) <- CCLECa$X
CCLECa <- CCLECa[,-1]
for (i in c(1 : 30)){
  a = 0
  for (j in c(1 : ncol(CCLECa))){
    if(CCLECa[i,j] != 0){
      a = a + 1
    }
  }
  Sum$CCLECa[i] = a
}

IhLine <- read.csv(file = "COSMIC_SIGNATURE_HNSC_local_cell_lines.csv")
rownames(IhLine) <- IhLine$X
IhLine <- IhLine[,-1]  
for (i in c(1:30)){
  a = 0
  for (j in c(1 : ncol(IhLine))){
    if(IhLine[i,j] != 0){
      a = a + 1
    }
  }
  Sum$IhLine[i] = a
}  

IhTu <- read.csv(file = "COSMIC_SIGNATURE_HNSC_local_tumors_20200120.csv")
rownames(IhTu) <- IhTu$X
IhTu <- IhTu[,-1]  
for (i in c(1:30)){
  a = 0
  for (j in c(1 : ncol(IhTu))){
    if(IhTu[i,j] != 0){
      a = a + 1
    }
  }
  Sum$IhTu[i] = a
} 

ICGC <- read.csv(file = "COSMIC_SIGNATURE_ICGC_tumors.csv")
rownames(ICGC) <- ICGC$X
ICGC <- ICGC[,-1]  
for (i in c(1:30)){
  a = 0
  for (j in c(1 : ncol(ICGC))){
    if(ICGC[i,j] != 0){
      a = a + 1
    }
  }
  Sum$ICGC[i] = a
} 

TCGA <- read_excel("COSMIC_SIGNATURE_TCGA-HNSC_tumors.xlsx")
TCGA <- TCGA[,-1]
TCGApatient <- read.csv(file = "data_bcr_clinical_data_patient.txt", sep = "\t")

TCGApatientAs <- TCGApatient %>% filter(Race.Category == "ASIAN")
AsPatientID <- as.vector(TCGApatientAs$Patient.Identifier)
AsPatientID <- AsPatientID[-2] #remove patient TCGA-BB-7872
TCGACOSMICAs <- TCGA %>% select(AsPatientID)

TCGApatientCa <- TCGApatient %>% filter(Race.Category == "WHITE")
CaPatientID1 <- as.vector(TCGApatientCa$Patient.Identifier)
CaPatientID2 <- as.vector(colnames(TCGA))
CaPatientID <- CaPatientID2[CaPatientID2 %in% CaPatientID1]
TCGACOSMICCa <- TCGA %>% select(CaPatientID)

for (i in c(1:30)){
  a = 0
  for (j in c(1 : ncol(TCGACOSMICAs))){
    if(TCGACOSMICAs[i,j] != 0){
      a = a + 1
    }
  }
  Sum$TCGAAs[i] = a
}

for (i in c(1:30)){
  a = 0
  for (j in c(1 : ncol(TCGACOSMICCa))){
    if(TCGACOSMICCa[i,j] != 0){
      a = a + 1
    }
  }
  Sum$TCGACa[i] = a
}

Sum$CCLEAs <- Sum$CCLEAs / 15  
Sum$CCLECa <- Sum$CCLECa /17
Sum$IhLine <- Sum$IhLine / 16
Sum$IhTu <- Sum$IhTu / 22
Sum$ICGC <- Sum$ICGC / 50
Sum$TCGAAs <- Sum$TCGAAs / 10
Sum$TCGACa <- Sum$TCGACa / 436

write.csv(Sum, file = "percentage_summary.csv")

