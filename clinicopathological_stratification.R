#Library preparation
library(TCGA2STAT)
library(tidyverse)
library(readxl)
library(expss)
library(gridExtra)

##Self defined function
#####
is.integer0 <- function(x){
  if (length(x) == 0L){
    return(0)
  }
  else {
    return(x)
  }
}
#####

#Site definition
#Site definition
#####
Oral_Cavity <- c("Alveolar Ridge", "Buccal Mucosa", "Floor of mouth", 
                 "Hard Palate", "Lip", "Oral Cavity", "Oral Tongue")
Oropharyngeal <- c("Oropharynx", "Tonsil", "Base of tongue")
Hypopharyngeal <- c("Hypopharynx")
Laryngeal <- c("Larynx")
site_definition <- site_definition %>% dplyr::select(Subsite, Major_site)
TCGA_Pat$MAJ_SITE <- vlookup(TCGA_Pat$PRIMARY_SITE, site_definition, "Major_site")
Inhouse_Pat <- read.csv("inhouse_patient_info.csv")
Inhouse_Pat$MAJ_SITE <- vlookup(Inhouse_Pat$PRIMARY_SITE, site_definition, "Major_site")
ICGC_Pat$MAJ_SITE <- "Oral Cavity"
#####

#Sites
#All tumors
#####
Sites <- c("Oral Cavity", "Oropharyngeal", "Laryngeal", "Hypopharyngeal", "Nasal cavity & Paranasal sinu", "Met")
count <- data.frame(row.names = Sites)
count$CaY <- sapply(Sites, function(x, y) {
  y = TCGA_Pat_Ca
  is.integer0(y %>%
    filter(MAJ_SITE == x) %>%
    summarise(n = length(unique(PATIENT_ID))) %>%
    .$n)
})
count$CaY[6] = 2
count$CaN <- 454 - count$CaY

#patient TCGA-KU-A6H7-06 TCGA-UF-A71A-06 have metastatic tumors
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-UF-A71A", "PRIMARY_SITE"]
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-UF-A71A", "HPV.status"]
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-KU-A6H7", "PRIMARY_SITE"]
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-KU-A6H7", "HPV.status"]

#count["Oral Cavity", "CaY"] <- count["Oral Cavity", "CaY"] + 1
#count["Oral Cavity", "CaN"] <- count["Oral Cavity", "CaN"] - 1
#count["Oropharyngeal", "CaY"] <- count["Oropharyngeal", "CaY"] + 1
#count["Oropharyngeal", "CaN"] <- count["Oropharyngeal", "CaN"] - 1

count$AsY <- sapply(Sites, function(x, y) {
  y = TCGA_Pat_As
  is.integer0(y %>%
                filter(MAJ_SITE == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$AsN <- 11 - count$AsY

count$IhY <- sapply(Sites, function(x, y) {
  y = Inhouse_Pat
  is.integer0(y %>%
                filter(MAJ_SITE == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IhY[6] = 1
count$IhN <- 22 - count$IhY
#count["Oral Cavity", "IhY"] <- count["Oral Cavity", "IhY"] + 1
#count["Oral Cavity", "IhN"] <- count["Oral Cavity", "IhN"] - 1


count$icY <- c(50, 0, 0, 0,0,0)
count$icN <- c(0, 50, 50, 50,50,50)

count$ComY <- count$AsY + count$IhY + count$icY
count$ComN <- count$AsN + count$IhN + count$icN

count$p1 <- 
  apply(count,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$p.value)
    }
  ) 
count$p2 <- apply(count,1,function(x) fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$p.value)

count$pctCaY <- round((count$CaY / 454) * 100, 2)
count$pctCaN <- round((count$CaN / 454) * 100, 2)
count$pctComY <- round((count$ComY / 83) * 100, 2)
count$pctComN <- round((count$ComN / 83) * 100, 2)

count$pctoutputCaY <- paste(count$CaY, " (", count$pctCaY, ")", sep = "")
count$pctoutputCaN <- paste(count$CaN, " (", count$pctCaN, ")", sep = "")
count$pctoutputComY <- paste(count$ComY, " (", count$pctComY, ")", sep = "")
count$pctoutputComN <- paste(count$ComN, " (", count$pctComN, ")", sep = "")

write_clip(count)
#####

#HPV-negative tumors
#####
Sites <- c("Oral Cavity", "Oropharyngeal", "Laryngeal", "Hypopharyngeal", "Nasal cavity & Paranasal sinu", "Met")
count <- data.frame(row.names = Sites)

count$CaY <- sapply(Sites, function(x, y) {
  y = TCGA_Pat_Ca_HPV_Neg
  is.integer0(y %>%
                filter(MAJ_SITE == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$CaY[6] = 1
count$CaN <- 352 - count$CaY

#patient TCGA-KU-A6H7-06 TCGA-UF-A71A-06 have metastatic tumors
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-UF-A71A", "PRIMARY_SITE"]
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-UF-A71A", "HPV.status"]
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-KU-A6H7", "PRIMARY_SITE"]
TCGA_Pat[TCGA_Pat$PATIENT_ID == "TCGA-KU-A6H7", "HPV.status"]

#count["Oral Cavity", "CaY"] <- count["Oral Cavity", "CaY"] + 1
#count["Oral Cavity", "CaN"] <- count["Oral Cavity", "CaN"] - 1

count$AsY <- sapply(Sites, function(x, y) {
  y = TCGA_Pat_As_HPV_Neg
  is.integer0(y %>%
                filter(MAJ_SITE == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$AsN <- 10 - count$AsY

count$IhY <- sapply(Sites, function(x, y) {
  y = Inhouse_Pat
  is.integer0(y %>%
                filter(MAJ_SITE == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IhY[6] = 1
count$IhN <- 22 - count$IhY
#count["Oral Cavity", "IhY"] <- count["Oral Cavity", "IhY"] + 1
#count["Oral Cavity", "IhN"] <- count["Oral Cavity", "IhN"] - 1


count$icY <- c(37, 0, 0, 0, 0, 0)
count$icN <- c(0, 37, 37, 37, 37, 37)

count$ComY <- count$AsY + count$IhY + count$icY
count$ComN <- count$AsN + count$IhN + count$icN

count$p1 <- 
  apply(count,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$p.value)
    }
  ) 
count$ODDS <- 
  apply(count,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$estimate)
    }
  ) 
count$pctCaY <- round((count$CaY / 352) * 100, 2)
count$pctCaN <- round((count$CaN / 352) * 100, 2)
count$pctComY <- round((count$ComY / 69) * 100, 2)
count$pctComN <- round((count$ComN / 69) * 100, 2)

count$pctoutputCaY <- paste(count$CaY, " (", count$pctCaY, ")", sep = "")
count$pctoutputCaN <- paste(count$CaN, " (", count$pctCaN, ")", sep = "")
count$pctoutputComY <- paste(count$ComY, " (", count$pctComY, ")", sep = "")
count$pctoutputComN <- paste(count$ComN, " (", count$pctComN, ")", sep = "")

write_clip(count)
#####



##Gender
#All tumor
#####
gender = c("Female", "Male")
count <- data.frame(row.names = gender)

count$CaY <- sapply(gender, function(x, y) {
  y = TCGA_Pat_Ca
  is.integer0(y %>%
                filter(SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$CaN <- 452 - count$CaY

count$AsY <- sapply(gender, function(x, y) {
  y = TCGA_Pat_As
  is.integer0(y %>%
                filter(SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$AsN <- 11 - count$AsY

count$IhY <- sapply(gender, function(x, y) {
  y = Inhouse_Pat
  is.integer0(y %>%
                filter(PATIENT_ID != "QM25" & SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IhN <- 21 - count$IhY

count$IcY <- sapply(gender, function(x, y) {
  y = ICGC_Pat
  is.integer0(y %>%
                filter(SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IcN <- 50 - count$IcY


count$ComY <- count$AsY + count$IhY + count$IcY
count$ComN <- count$AsN + count$IhN + count$IcN

count$p1 <- 
  apply(count,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$p.value)
    }
  ) 
count$p2 <- apply(count,1,function(x) fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$p.value)

count$pctCaY <- round((count$CaY / 452) * 100, 2)
count$pctCaN <- round((count$CaN / 452) * 100, 2)
count$pctComY <- round((count$ComY / 82) * 100, 2)
count$pctComN <- round((count$ComN / 82) * 100, 2)

count$pctoutputCaY <- paste(count$CaY, " (", count$pctCaY, ")", sep = "")
count$pctoutputCaN <- paste(count$CaN, " (", count$pctCaN, ")", sep = "")
count$pctoutputComY <- paste(count$ComY, " (", count$pctComY, ")", sep = "")
count$pctoutputComN <- paste(count$ComN, " (", count$pctComN, ")", sep = "")

write_clip(count)
#####

#HPV-negative tumor
#####
gender = c("Female", "Male")
count <- data.frame(row.names = gender)

count$CaY <- sapply(gender, function(x, y) {
  y = TCGA_Pat_Ca
  is.integer0(y %>%
                filter(HPV.status == "Neg" & SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$CaN <- 351 - count$CaY

count$AsY <- sapply(gender, function(x, y) {
  y = TCGA_Pat_As_HPV_Neg
  is.integer0(y %>%
                filter(SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$AsN <- 10 - count$AsY

count$IhY <- sapply(gender, function(x, y) {
  y = Inhouse_Pat
  is.integer0(y %>%
                filter(PATIENT_ID != "QM25" & SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IhN <- 21 - count$IhY

count$IcY <- sapply(gender, function(x, y) {
  y = ICGC_Pat
  is.integer0(y %>%
                filter(HPV.status == "Neg", SEX == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IcN <- 37 - count$IcY


count$ComY <- count$AsY + count$IhY + count$IcY
count$ComN <- count$AsN + count$IhN + count$IcN

count$p1 <- 
  apply(count,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$p.value)
    }
  ) 
count$p2 <- apply(count,1,function(x) fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$p.value)

count$pctCaY <- round((count$CaY / 351) * 100, 2)
count$pctCaN <- round((count$CaN / 351) * 100, 2)
count$pctComY <- round((count$ComY / 68) * 100, 2)
count$pctComN <- round((count$ComN / 68) * 100, 2)

count$pctoutputCaY <- paste(count$CaY, " (", count$pctCaY, ")", sep = "")
count$pctoutputCaN <- paste(count$CaN, " (", count$pctCaN, ")", sep = "")
count$pctoutputComY <- paste(count$ComY, " (", count$pctComY, ")", sep = "")
count$pctoutputComN <- paste(count$ComN, " (", count$pctComN, ")", sep = "")

write_clip(count)
#####



##PN
#All tumor
#####
PN = c("YES", "NO")
count <- data.frame(row.names = PN)

count$CaY <- sapply(PN, function(x, y) {
  y = TCGA_Pat_Ca
  is.integer0(y %>%
                filter(PERINEURAL_INVASION == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$CaN <- 312 - count$CaY

count$AsY <- sapply(PN, function(x, y) {
  y = TCGA_Pat_As
  is.integer0(y %>%
                filter(PERINEURAL_INVASION == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$AsN <- 9 - count$AsY

count$IhY <- sapply(PN, function(x, y) {
  y = Inhouse_Pat
  is.integer0(y %>%
                filter(PERINEURAL_INVASION == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IhN <- 22 - count$IhY



count$ComY <- count$AsY + count$IhY 
count$ComN <- count$AsN + count$IhN 

count$p1 <- 
  apply(count,1,function(x)
    if (fisher.test(matrix(x[c(1,2,7,8)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,7,8)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,7,8)],nrow=2), alternative = "greater")$p.value)
    }
  ) 
count$p2 <- apply(count,1,function(x) fisher.test(matrix(x[c(1,2,7,8)],nrow=2))$p.value)

count$pctCaY <- round((count$CaY / 312) * 100, 2)
count$pctCaN <- round((count$CaN / 312) * 100, 2)
count$pctComY <- round((count$ComY / 31) * 100, 2)
count$pctComN <- round((count$ComN / 31) * 100, 2)

count$pctoutputCaY <- paste(count$CaY, " (", count$pctCaY, ")", sep = "")
count$pctoutputCaN <- paste(count$CaN, " (", count$pctCaN, ")", sep = "")
count$pctoutputComY <- paste(count$ComY, " (", count$pctComY, ")", sep = "")
count$pctoutputComN <- paste(count$ComN, " (", count$pctComN, ")", sep = "")

write_clip(count)
#####

#HPV-negative tumor
#####
PN = c("YES", "NO")
count <- data.frame(row.names = PN)

count$CaY <- sapply(PN, function(x, y) {
  y = TCGA_Pat_Ca_HPV_Neg
  is.integer0(y %>%
                filter(PERINEURAL_INVASION == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$CaN <- 250 - count$CaY

count$AsY <- sapply(PN, function(x, y) {
  y = TCGA_Pat_As_HPV_Neg
  is.integer0(y %>%
                filter(PERINEURAL_INVASION == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$AsN <- 8 - count$AsY

count$IhY <- sapply(PN, function(x, y) {
  y = Inhouse_Pat
  is.integer0(y %>%
                filter(PERINEURAL_INVASION == x) %>%
                summarise(n = length(unique(PATIENT_ID))) %>%
                .$n)
})
count$IhN <- 22 - count$IhY


count$ComY <- count$AsY + count$IhY
count$ComN <- count$AsN + count$IhN

count$p1 <- 
  apply(count,1,function(x)
    if (fisher.test(matrix(x[c(1,2,7,8)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,7,8)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,7,8)],nrow=2), alternative = "greater")$p.value)
    }
  ) 
count$p2 <- apply(count,1,function(x) fisher.test(matrix(x[c(1,2,7,8)],nrow=2))$p.value)

count$pctCaY <- round((count$CaY / 250) * 100, 2)
count$pctCaN <- round((count$CaN / 250) * 100, 2)
count$pctComY <- round((count$ComY / 30) * 100, 2)
count$pctComN <- round((count$ComN / 30) * 100, 2)

count$pctoutputCaY <- paste(count$CaY, " (", count$pctCaY, ")", sep = "")
count$pctoutputCaN <- paste(count$CaN, " (", count$pctCaN, ")", sep = "")
count$pctoutputComY <- paste(count$ComY, " (", count$pctComY, ")", sep = "")
count$pctoutputComN <- paste(count$ComN, " (", count$pctComN, ")", sep = "")

write_clip(count)
#####






























#####























df1 <- data.frame()
Sites <- unique(c(as.character(TCGA_Pat$MAJ_SITE), as.character(Inhouse_Pat$MAJ_SITE)))
Sites <- Sites[-14]
for (i in c(1:13)){
  site = Sites[i]
  pl1 <- TCGA_Pat_Ca %>% filter(MAJ_SITE == site) %>% .$PATIENT_ID
  temp1 <- TCGA_X_Ca %>% 
    filter(Patient_ID %in% pl1) %>% 
    filter(Hugo_Symbol %in% cancer_gene_census_combined) %>%
    group_by(Hugo_Symbol) %>%
    summarise(n = length(unique(Tumor_Sample_Barcode)), Site = site, Cohort = "TCGA Caucasian")
  pl2 <- TCGA_Pat_As %>% filter(MAJ_SITE == site) %>% .$PATIENT_ID
  temp2 <- TCGA_X_As %>% 
    filter(Patient_ID %in% pl2) %>% 
    filter(Hugo_Symbol %in% cancer_gene_census_combined) %>%
    group_by(Hugo_Symbol) %>%
    summarise(n = length(unique(Tumor_Sample_Barcode)), Site = site, Cohort = "TCGA Asian")
  pl3 <- Inhouse_Pat %>% filter(MAJ_SITE == site) %>% .$PATIENT_ID
  temp3 <- Inhouse_tumor_X %>% 
    filter(Tumor_Sample_Barcode %in% pl3) %>% 
    filter(Hugo_Symbol %in% cancer_gene_census_combined) %>%
    group_by(Hugo_Symbol) %>%
    summarise(n = length(unique(Tumor_Sample_Barcode)), Site = site, Cohort = "Inhouse")
  df1 <- rbind(df, temp1, temp2, temp3)
}
#####

#Construct fisher matrix
#####
output = data.frame()
#all mutations
i = 11
for (i in c(1:12)) {
  site = Sites[i]
  FisherMat <- data.frame(row.names = cancer_gene_census_combined)
  
  FisherMat$CauMut <- sapply(cancer_gene_census_combined, function(x) 
    is.integer0(df1 %>% filter(Site == site & Cohort == "TCGA Caucasian" & Hugo_Symbol == x) %>% .$n))
  FisherMat$CauUnmut <- is.integer0(count %>% filter(MAJ_SITE == site & Cohort == "TCGA Caucasian") %>% .$n) - FisherMat$CauMut
  
  FisherMat$AsMut <- sapply(cancer_gene_census_combined, function(x) 
    is.integer0(df1 %>% filter(Site == site & Cohort == "TCGA Asian" & Hugo_Symbol == x) %>% .$n))
  FisherMat$AsUnmut <- is.integer0(count %>% filter(MAJ_SITE == site & Cohort == "TCGA Asian") %>% .$n) - FisherMat$AsMut
  
  FisherMat$IhMut <- sapply(cancer_gene_census_combined, function(x) 
    is.integer0(df1 %>% filter(Site == site & Cohort == "Inhouse" & Hugo_Symbol == x) %>% .$n))
  FisherMat$IhUnmut <- is.integer0(count %>% filter(MAJ_SITE == site & Cohort == "Inhouse") %>% .$n) - FisherMat$IhMut
  
  #combine two
  FisherMat$TcIhMut <- FisherMat$AsMut + FisherMat$IhMut
  FisherMat$TcIhUnmut <- FisherMat$AsUnmut + FisherMat$IhUnmut
  
  FisherMat$p2 <- apply(FisherMat,1,function(x) fisher.test(matrix(x[c(1,2,7,8)],nrow=2))$p.value)
  FisherMat$p1 <- 
    apply(FisherMat,1,function(x)
      if (fisher.test(matrix(x[c(1,2,7,8)],nrow=2))$estimate < 1){
        return(fisher.test(matrix(x[c(1,2,7,8)],nrow=2), alternative = "less")$p.value)
      }
      else {
        return(fisher.test(matrix(x[c(1,2,7,8)],nrow=2), alternative = "greater")$p.value)
      }
    )
  FisherMat$Gene <- rownames(FisherMat)
  FisherMat$Site <- site
  l1 <- FisherMat %>% filter(p1 < 0.05 | p1 == 0.05)
  l2 <- FisherMat %>% filter(p2 < 0.05 | p2 == 0.05)
  #l1$Site <- site
  #l2$Site <- site
  output <- rbind(output, l1, l2)
}

write_clip(output)
#####

#Gender
#####
gen_ca <- TCGA_Pat_Ca %>%
  group_by(SEX) %>%
  summarise(Ca = length(unique(PATIENT_ID)))

gen_ca_hpv_neg <- TCGA_Pat_Ca_HPV_Neg %>%
  group_by(SEX) %>%
  summarise(Ca_hpv_neg = length(unique(PATIENT_ID)))

gen_as <- TCGA_Pat_As %>%
  group_by(SEX) %>%
  summarise(TCGA_as= length(unique(PATIENT_ID)))

gen_ih <- Inhouse_Pat %>%
  group_by(SEX) %>%
  summarise(Ih = length(unique(PATIENT_ID)))

gen_ic <- ICGC_Pat %>%
  group_by(SEX) %>%
  summarise(ICGC = length(unique(PATIENT_ID)))

gen_ic_hpv_neg <- ICGC_Pat %>%
  filter(HPV.status == "Neg") %>%
  group_by(SEX) %>%
  summarise(ICGC_HPV_Neg = length(unique(PATIENT_ID)))

count <- cbind(gen_ca, gen_as, gen_ca_hpv_neg, gen_ic, gen_ic_hpv_neg, gen_ih)
write_clip(count)
#####

#Perineural invasion
#####
pn_ca <- TCGA_Pat_Ca %>%
  group_by(PERINEURAL_INVASION) %>%
  filter(PERINEURAL_INVASION != "[Not Available]") %>%
  summarise(Ca = length(unique(PATIENT_ID)))

pn_ca_hpv_neg <- TCGA_Pat_Ca_HPV_Neg %>%
  group_by(PERINEURAL_INVASION) %>%
  filter(PERINEURAL_INVASION != "[Not Available]") %>%
  summarise(Ca_hpv_neg = length(unique(PATIENT_ID)))

pn_as <- TCGA_Pat_As %>%
  group_by(PERINEURAL_INVASION) %>%
  filter(PERINEURAL_INVASION != "[Not Available]") %>%
  summarise(TCGA_as= length(unique(PATIENT_ID)))

pn_ih <- Inhouse_Pat %>%
  group_by(PERINEURAL_INVASION) %>%
  summarise(Ih = length(unique(PATIENT_ID)))

count <- cbind(pn_ca, pn_as, pn_ca_hpv_neg, pn_ih)
write_clip(count)
#####

#Age
