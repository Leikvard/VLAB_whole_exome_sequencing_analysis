#Library preparation
library(tidyverse)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(expss)


##Self defined function
##########################################################################################
is.integer0 <- function(x){
  if (length(x) == 0L){
    return(0)
  }
  else {
    return(x)
  }
}

#---------------------------------------------------------------------------------------------

#nonsyn mutations
#HPV-negative
#Oral Cavity
#DMG
##########################################################################################
#Caucasian
pl_ca <- TCGA_Pat_Ca_HPV_Neg %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_ca <- as.character(pl_ca)
CauCount <- TCGA_X_Ca_HPV_Neg %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode != "TCGA-UF-A71A-06") %>%
  filter(Patient_ID %in% pl_ca) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
CauCount <- as.data.frame(CauCount)

#Asian
pl_as <- TCGA_Pat_As_HPV_Neg %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_as <- as.character(pl_as)
AsCount <- TCGA_X_As_HPV_Neg %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Patient_ID %in% pl_as) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
AsCount <- as.data.frame(AsCount)

#icgc hpv-neg ns
icgcCount <- ICGC_HPV_Neg %>% 
  filter(Variant_Classification %in% NS) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
icgcCount <- as.data.frame(icgcCount)

#inhouse
pl_ih <- Inhouse_Pat %>% filter(MAJ_SITE == "Oral Cavity") %>% .$PATIENT_ID
pl_ih <- as.character(pl_ih)
IhCount <- Inhouse_tumor_X %>% 
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode %in% pl_ih) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
IhCount <- as.data.frame(IhCount)

#Construct fisher matrix
g <- unique(as.character(TCGA_X$Hugo_Symbol))

FisherMat <- data.frame(row.names = g)

#mutations
FisherMat$CauMut <- sapply(g, function(x) is.integer0(CauCount[CauCount$Hugo_Symbol == x,2]))
FisherMat$CauUnmut <- length(unique(TCGA_X_Ca_HPV_Neg %>% filter(Patient_ID %in% pl_ca & Tumor_Sample_Barcode != "TCGA-UF-A71A-06") %>% .$Tumor_Sample_Barcode)) - FisherMat$CauMut

FisherMat$TAsMut <- sapply(g, function(x) is.integer0(AsCount[AsCount$Hugo_Symbol == x,2]))
FisherMat$TAsUnmut <- length(unique(TCGA_X_As_HPV_Neg %>%
                                      filter(Patient_ID %in% pl_as) %>% .$Tumor_Sample_Barcode)) - FisherMat$TAsMut

FisherMat$InhouseMut <- sapply(g, function(x) is.integer0(IhCount[IhCount$Hugo_Symbol == x,2]))
FisherMat$InhouseUnmut <- length(unique(Inhouse_tumor_X %>%
                                          filter(Tumor_Sample_Barcode %in% pl_ih) %>% .$Tumor_Sample_Barcode)) - FisherMat$InhouseMut

FisherMat$icgcMut <- sapply(g, function(x) is.integer0(icgcCount[icgcCount$Hugo_Symbol == x,2]))
FisherMat$icgcUnmut <- n_ICGC_HPV_Neg - FisherMat$icgcMut

#combine two
FisherMat$TcIhMut <- FisherMat$TAsMut + FisherMat$InhouseMut
FisherMat$TcIhUnmut <- FisherMat$TAsUnmut + FisherMat$InhouseUnmut

#combine three
FisherMat$TcIhIcMut <- FisherMat$TAsMut + FisherMat$InhouseMut + FisherMat$icgcMut
FisherMat$TcIhIcUnmut <- FisherMat$TAsUnmut + FisherMat$InhouseUnmut + FisherMat$icgcUnmut

#calculate p_value
#FisherMat$p2_2As <- apply(FisherMat,1,function(x) fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$p.value)
#FisherMat$p2_3As <- apply(FisherMat,1,function(x) fisher.test(matrix(x[c(1,2,11,12)],nrow=2))$p.value)
FisherMat$p1_2As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$p.value)
    }
  )
FisherMat$p1_3As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,11,12)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "greater")$p.value)
    }
  )
FisherMat$odds_2As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$estimate)
    }
  )
FisherMat$odds_3As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,11,12)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "greater")$estimate)
    }
  )
FisherMat$Gene <- rownames(FisherMat)
write_clip(FisherMat %>% filter(p1_2As < 0.05 | p1_2As == 0.05))
write_clip(FisherMat %>% filter(p1_3As < 0.05 | p1_3As == 0.05))


#two asian
dmg2 <- FisherMat %>% filter(p1_2As < 0.05 | p1_2As == 0.05) %>% .$Gene

gene.df <- bitr(dmg2, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
#write.csv(ego, "test.csv")
ego_sig <- ego[ego$pvalue < 0.05, asis =T]
barplot(ego_sig,x='Count', color = 'pvalue', showCategory = 30)

kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") # setReadable to see the matched genes
#head(kk2)
#write.csv(kk2, "common_kegg_output.csv")
kk2_sig <- kk2[kk2$pvalue < 0.05, asis =T]
barplot(kk2_sig,x='Count', color = 'pvalue', showCategory = 20)

#three asian
dmg3 <- FisherMat %>% filter(p1_3As < 0.05 | p1_3As == 0.05) %>% .$Gene
gene.df <- bitr(dmg3, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
#write.csv(ego, "test.csv")
ego_sig <- ego[ego$pvalue < 0.05, asis =T]
barplot(ego_sig,x='Count', color = 'pvalue', showCategory = 30)

kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") # setReadable to see the matched genes
#head(kk2)

#write.csv(kk2, "common_kegg_output.csv")
kk2_sig <- kk2[kk2$pvalue < 0.05, asis =T]
barplot(kk2_sig,x='Count', color = 'pvalue', showCategory = 20)
#####


#nonsyn mutations
#HPV-positive
#Oral Cavity
#DMG
##########################################################################################
#Caucasian hpv-pos
pl_ca <- TCGA_Pat_Ca_HPV_Pos %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_ca <- as.character(pl_ca)
CauCount <- TCGA_X_Ca_HPV_Pos %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Patient_ID %in% pl_ca) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
CauCount <- as.data.frame(CauCount)

#icgc hpv-pos ns
icgcCount <- ICGC_HPV_Pos %>% 
  filter(Variant_Classification %in% NS) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
icgcCount <- as.data.frame(icgcCount)

#Construct fisher matrix
g <- unique(as.character(TCGA_X$Hugo_Symbol))

FisherMat <- data.frame(row.names = g)

#mutations
FisherMat$CauMut <- sapply(g, function(x) is.integer0(CauCount[CauCount$Hugo_Symbol == x,2]))
FisherMat$CauUnmut <- length(unique(TCGA_X_Ca_HPV_Pos %>% filter(Patient_ID %in% pl_ca) %>% .$Tumor_Sample_Barcode)) - FisherMat$CauMut

FisherMat$icgcMut <- sapply(g, function(x) is.integer0(icgcCount[icgcCount$Hugo_Symbol == x,2]))
FisherMat$icgcUnmut <- n_ICGC_HPV_Pos - FisherMat$icgcMut

#calculate p_value
FisherMat$p1 <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1:4)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "greater")$p.value)
    }
  )

FisherMat$odds <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1:4)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "greater")$estimate)
    }
  )

FisherMat$Gene <- rownames(FisherMat)

#Map to pathways
dmg <- FisherMat %>% filter(p1 < 0.05 | p1 == 0.05) %>% .$Gene

gene.df <- bitr(dmg, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
#write.csv(ego, "test.csv")
ego_sig <- ego[ego$pvalue < 0.05, asis =T]
barplot(ego_sig,x='Count', color = 'pvalue', showCategory = 30)

kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") # setReadable to see the matched genes
#head(kk2)
#write.csv(kk2, "common_kegg_output.csv")
kk2_sig <- kk2[kk2$pvalue < 0.05, asis =T]
barplot(kk2_sig,x='Count', color = 'pvalue', showCategory = 20)
#####

#convert gene names to hugo symbol
#####
#convert tcga
x <- read.csv("./Environment/tcga_hugo.csv", header = FALSE)
x <- as.data.frame(x)
x <- x[-c(1,2),]
x <- x %>% dplyr::select(c(1:3))
colnames(x) <- c("Input", "Match_type", "Hugo")
x$Input <- as.character(x$Input)
x$Hugo <- as.character(x$Hugo)
y <- x %>% filter(Input == Hugo)
y <- y[!duplicated(y$Input),]
y$Hugo[duplicated(y$Hugo)]
x <- x %>% filter(Input != Hugo) %>% filter(!Input %in% y$Input) %>% filter(Hugo != "")
x$Input[duplicated(x$Input)]
x <- x[-which(x$Input == "NOV" & x$Hugo != "CCN3"),] #remove NOV that is translated to RPL10 or PLXNA1, because TCGA_X NOV is on chr8, should be CCN3
x <- x[-which(x$Input == "QARS" & x$Hugo != "QARS1"),] #QARS should be QARS1
x <- x[-which(x$Input == "MUM1" & x$Hugo != "PWWP3A"),]
x <- x[-which(x$Input == "C11orf48" & x$Hugo != "LBHD1"),]
x <- x[-which(x$Input == "CSRP2BP" & x$Hugo != "KAT14"),]
x <- x[-which(x$Input == "GIF" & x$Hugo != "CBLIF"),]
x <- x[-which(x$Input == "B3GNT1" & x$Hugo != "B4GAT1"),]
x <- x[-which(x$Input == "DEC1"),] #DEC1 is a valid Hugo symbol
x <- x[-which(x$Input == "CXXC11" & x$Hugo != "RTP5"),]
x <- x[-which(x$Input == "SEPT2" & x$Hugo != "SEPTIN2"),]
x <- x[-which(x$Input == "SARS" & x$Hugo != "SARS1"),]
x <- x[-which(x$Input == "AGPAT9" & x$Hugo != "GPAT3"),]
x <- x[-which(x$Input == "LOR"),] #LOR is a valid hugo symbol
x <- x[-which(x$Input == "ISPD")[1],]
x$Hugo[duplicated(x$Hugo)]
dic <- rbind(x,y)
dic$Input[duplicated(dic$Input)]
dic$Hugo[duplicated(dic$Hugo)]
#all checked, inputs with same hugo are actually same genes
TCGA_X_copy <- TCGA_X
TCGA_X_copy$Hugo_Symbol <- sapply(TCGA_X_copy$Hugo_Symbol, function(x) {
  if (x %in% dic$Input) {
    return(vlookup(x, dic, "Hugo"))
  }
  else {
    return(x)
  }
})
length(unique(TCGA_X_copy$Hugo_Symbol))
TCGA_X_Ca_HPV_Pos_copy <- TCGA_X_copy %>% filter(Patient_ID %in% TCGA_Pat_Ca_HPV_Pos$PATIENT_ID)
TCGA_X_Ca_HPV_Neg_copy <- TCGA_X_copy %>% filter(Patient_ID %in% TCGA_Pat_Ca_HPV_Neg$PATIENT_ID)
TCGA_X_As_HPV_Neg_copy <- TCGA_X_copy %>% filter(Patient_ID %in% TCGA_Pat_As_HPV_Neg$PATIENT_ID)


#convert inhouse
x <- read.csv("./Environment/Inhouse_hugo.csv", header = FALSE)
x <- as.data.frame(x)
x <- x[-c(1,2),]
x <- x %>% dplyr::select(c(1:3))
colnames(x) <- c("Input", "Match_type", "Hugo")
x$Input <- as.character(x$Input)
x$Hugo <- as.character(x$Hugo)
y <- x %>% filter(Input == Hugo)
y <- y[!duplicated(y$Input),]
x <- x %>% filter(Input != Hugo) %>% filter(!Input %in% y$Input) %>% filter(Hugo != "")
x$Input[duplicated(x$Input)]
  #MUM1 has two approved symbol, since MUM1 in Inhouse MAF file is on chr19, so it should be PWWP3A
x <- x[-95,] #remove MUM1 that is translated to IRF4
x$Hugo[duplicated(x$Hugo)]
  #CHDC2 and CXorf22 are both converted to CFAP47, they are all on X, position close to each other, probably same gene
y$Hugo[duplicated(y$Hugo)]
dic <- rbind(x,y)
dic$Input[duplicated(dic$Input)]
dic$Hugo[duplicated(dic$Hugo)]
  #all checked, inputs with same hugo are actually same genes
Inhouse_copy <- Inhouse_tumor_X
#Inhouse_copy$Hugo_Symbol <- vlookup(Inhouse_copy$Hugo_Symbol, dic, "Hugo")
Inhouse_copy$Hugo_Symbol <- sapply(Inhouse_copy$Hugo_Symbol, function(x) {
  if (x %in% dic$Input) {
    return(vlookup(x, dic, "Hugo"))
  }
  else {
    return(x)
  }
})

#convert ICGC
x <- read.csv("./Environment/icgc_hugo.csv", header = FALSE)
x <- as.data.frame(x)
x <- x[-c(1,2),]
x <- x %>% dplyr::select(c(1:3))
colnames(x) <- c("Input", "Match_type", "Hugo")
x$Input <- as.character(x$Input)
x$Hugo <- as.character(x$Hugo)
y <- x %>% filter(Input == Hugo)
y <- y[!duplicated(y$Input),]
x <- x %>% filter(Input != Hugo) %>% filter(!Input %in% y$Input) %>% filter(Hugo != "")
x$Input[duplicated(x$Input)]
#MUM1 has two approved symbol, since MUM1 in Inhouse MAF file is on chr19, so it should be PWWP3A
x <- x[-which(x$Input == "MLL2" & x$Hugo != "KMT2D"),]
x <- x[-which(x$Input == "MLL4" & x$Hugo != "KMT2B"),]
x <- x[-which(x$Input == "ODZ3" & x$Hugo != "TENM3"),]
x$Hugo[duplicated(x$Hugo)]
#CHDC2 and CXorf22 are both converted to CFAP47, they are all on X, position close to each other, probably same gene
y$Hugo[duplicated(y$Hugo)]
dic <- rbind(x,y)
dic$Input[duplicated(dic$Input)]
dic$Hugo[duplicated(dic$Hugo)]
#all checked, inputs with same hugo are actually same genes
ICGC_copy <- ICGC
ICGC_copy$Hugo_Symbol <- sapply(ICGC_copy$Hugo_Symbol, function(x) {
  if (x %in% dic$Input) {
    return(vlookup(x, dic, "Hugo"))
  }
  else {
    return(x)
  }
})
length(unique(ICGC_copy$Hugo_Symbol))
sum(is.na(ICGC_copy$Hugo_Symbol))
ICGC_HPV_Pos_copy <- ICGC_copy %>%
  filter(`Tumor_Sample_Barcode` %in% ICGC_Pat[ICGC_Pat$HPV.status == "Pos",]$PATIENT_ID)
ICGC_HPV_Neg_copy <- ICGC_copy %>%
  filter(`Tumor_Sample_Barcode` %in% ICGC_Pat[ICGC_Pat$HPV.status == "Neg",]$PATIENT_ID)

#####

#nonsyn mutations
#Hugo symbol converted
#HPV-negative
#Oral Cavity
#DMG
##########################################################################################
#Caucasian
pl_ca <- TCGA_Pat_Ca_HPV_Neg %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_ca <- as.character(pl_ca)
CauCount <- TCGA_X_Ca_HPV_Neg_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode != "TCGA-UF-A71A-06") %>%
  filter(Patient_ID %in% pl_ca) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
CauCount <- as.data.frame(CauCount)

#Asian
pl_as <- TCGA_Pat_As_HPV_Neg %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_as <- as.character(pl_as)
AsCount <- TCGA_X_As_HPV_Neg_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Patient_ID %in% pl_as) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
AsCount <- as.data.frame(AsCount)

#icgc hpv-neg ns
icgcCount <- ICGC_HPV_Neg_copy %>% 
  filter(Variant_Classification %in% NS) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
icgcCount <- as.data.frame(icgcCount)

#inhouse
pl_ih <- Inhouse_Pat %>% filter(MAJ_SITE == "Oral Cavity") %>% .$PATIENT_ID
pl_ih <- as.character(pl_ih)
IhCount <- Inhouse_copy %>% 
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode %in% pl_ih) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
IhCount <- as.data.frame(IhCount)

#Construct fisher matrix
g <- unique(c(as.character(TCGA_X_copy$Hugo_Symbol), as.character(Inhouse_copy$Hugo_Symbol), as.character(ICGC_copy$Hugo_Symbol)))

FisherMat <- data.frame(row.names = g)

#mutations
FisherMat$CauMut <- sapply(g, function(x) is.integer0(CauCount[CauCount$Hugo_Symbol == x,2]))
FisherMat$CauUnmut <- length(unique(TCGA_X_Ca_HPV_Neg %>% filter(Patient_ID %in% pl_ca & Tumor_Sample_Barcode != "TCGA-UF-A71A-06") %>% .$Tumor_Sample_Barcode)) - FisherMat$CauMut

FisherMat$TAsMut <- sapply(g, function(x) is.integer0(AsCount[AsCount$Hugo_Symbol == x,2]))
FisherMat$TAsUnmut <- length(unique(TCGA_X_As_HPV_Neg %>%
                                      filter(Patient_ID %in% pl_as) %>% .$Tumor_Sample_Barcode)) - FisherMat$TAsMut

FisherMat$InhouseMut <- sapply(g, function(x) is.integer0(IhCount[IhCount$Hugo_Symbol == x,2]))
FisherMat$InhouseUnmut <- length(unique(Inhouse_tumor_X %>%
                                          filter(Tumor_Sample_Barcode %in% pl_ih) %>% .$Tumor_Sample_Barcode)) - FisherMat$InhouseMut

FisherMat$icgcMut <- sapply(g, function(x) is.integer0(icgcCount[icgcCount$Hugo_Symbol == x,2]))
FisherMat$icgcUnmut <- n_ICGC_HPV_Neg - FisherMat$icgcMut

#combine two
FisherMat$TcIhMut <- FisherMat$TAsMut + FisherMat$InhouseMut
FisherMat$TcIhUnmut <- FisherMat$TAsUnmut + FisherMat$InhouseUnmut

#combine three
FisherMat$TcIhIcMut <- FisherMat$TAsMut + FisherMat$InhouseMut + FisherMat$icgcMut
FisherMat$TcIhIcUnmut <- FisherMat$TAsUnmut + FisherMat$InhouseUnmut + FisherMat$icgcUnmut

#calculate p_value
#FisherMat$p2_2As <- apply(FisherMat,1,function(x) fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$p.value)
#FisherMat$p2_3As <- apply(FisherMat,1,function(x) fisher.test(matrix(x[c(1,2,11,12)],nrow=2))$p.value)
FisherMat$p1_2As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$p.value)
    }
  )
FisherMat$p1_3As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,11,12)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "greater")$p.value)
    }
  )
FisherMat$odds_2As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,9,10)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,9,10)],nrow=2), alternative = "greater")$estimate)
    }
  )
FisherMat$odds_3As <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1,2,11,12)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1,2,11,12)],nrow=2), alternative = "greater")$estimate)
    }
  )
FisherMat$Gene <- rownames(FisherMat)
write_clip(FisherMat %>% filter(p1_2As < 0.05 | p1_2As == 0.05))
write_clip(FisherMat %>% filter(p1_3As < 0.05 | p1_3As == 0.05))


#two asian
dmg2 <- FisherMat %>% filter(p1_2As < 0.05 | p1_2As == 0.05) %>% .$Gene

gene.df <- bitr(dmg2, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
#write.csv(ego, "test.csv")
ego_sig <- ego[ego$pvalue < 0.05, asis =T]
barplot(ego_sig,x='Count', color = 'pvalue', showCategory = 30)

kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") # setReadable to see the matched genes
#head(kk2)
#write.csv(kk2, "common_kegg_output.csv")
kk2_sig <- kk2[kk2$pvalue < 0.05, asis =T]
barplot(kk2_sig,x='Count', color = 'pvalue', showCategory = 20)

#three asian
dmg3 <- FisherMat %>% filter(p1_3As < 0.05 | p1_3As == 0.05) %>% .$Gene
gene.df <- bitr(dmg3, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
#write.csv(ego, "test.csv")
ego_sig <- ego[ego$pvalue < 0.05, asis =T]
barplot(ego_sig,x='Count', color = 'pvalue', showCategory = 30)

kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") # setReadable to see the matched genes
#head(kk2)

#write.csv(kk2, "common_kegg_output.csv")
kk2_sig <- kk2[kk2$pvalue < 0.05, asis =T]
barplot(kk2_sig,x='Count', color = 'pvalue', showCategory = 20)
#####

#nonsyn mutations
#Hugo symbol converted
#HPV-positive
#Oral Cavity
#DMG
##########################################################################################
#Caucasian hpv-pos
pl_ca <- TCGA_Pat_Ca_HPV_Pos %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_ca <- as.character(pl_ca)
CauCount <- TCGA_X_Ca_HPV_Pos_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Patient_ID %in% pl_ca) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
CauCount <- as.data.frame(CauCount)

#icgc hpv-pos ns
icgcCount <- ICGC_HPV_Pos_copy %>% 
  filter(Variant_Classification %in% NS) %>%
  group_by(Hugo_Symbol) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode)))
icgcCount <- as.data.frame(icgcCount)

#Construct fisher matrix
g <- unique(c(as.character(TCGA_X_copy$Hugo_Symbol), as.character(ICGC_copy$Hugo_Symbol)))

FisherMat <- data.frame(row.names = g)

#mutations
FisherMat$CauMut <- sapply(g, function(x) is.integer0(CauCount[CauCount$Hugo_Symbol == x,2]))
FisherMat$CauUnmut <- length(unique(TCGA_X_Ca_HPV_Pos %>% filter(Patient_ID %in% pl_ca) %>% .$Tumor_Sample_Barcode)) - FisherMat$CauMut

FisherMat$icgcMut <- sapply(g, function(x) is.integer0(icgcCount[icgcCount$Hugo_Symbol == x,2]))
FisherMat$icgcUnmut <- n_ICGC_HPV_Pos - FisherMat$icgcMut

#calculate p_value
FisherMat$p1 <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1:4)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "greater")$p.value)
    }
  )

FisherMat$odds <- 
  apply(FisherMat,1,function(x)
    if (fisher.test(matrix(x[c(1:4)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "greater")$estimate)
    }
  )

FisherMat$Gene <- rownames(FisherMat)
write_clip(FisherMat %>% filter(p1 < 0.05 | p1 == 0.05))
#Map to pathways
dmg <- FisherMat %>% filter(p1 < 0.05 | p1 == 0.05) %>% .$Gene

gene.df <- bitr(dmg, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
#write.csv(ego, "test.csv")
ego_sig <- ego[ego$pvalue < 0.05, asis =T]
barplot(ego_sig,x='Count', color = 'pvalue', showCategory = 30)

kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") # setReadable to see the matched genes
#head(kk2)
#write.csv(kk2, "common_kegg_output.csv")
kk2_sig <- kk2[kk2$pvalue < 0.05, asis =T]
barplot(kk2_sig,x='Count', color = 'pvalue', showCategory = 20)
#####

