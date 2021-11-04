#Library preparation
library(tidyverse)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(expss)

#all cds feature table from hg38 assembly, downloaded from ncbi
cds <- data.table::fread("./Environment/GCF_000001405.39_GRCh38.p13_feature_table.txt", sep = "\t")
colnames(cds)[1] <- "feature"

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
#Hugo symbol converted
#HPV-negative
#Oral Cavity
#DMG_fisher analysis
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
FisherMat$CauUnmut <- length(unique(TCGA_X_Ca_HPV_Neg %>% 
                                      filter(Patient_ID %in% pl_ca & Tumor_Sample_Barcode != "TCGA-UF-A71A-06") %>% .$Tumor_Sample_Barcode)) - FisherMat$CauMut

FisherMat$TAsMut <- sapply(g, function(x) is.integer0(AsCount[AsCount$Hugo_Symbol == x,2]))
FisherMat$TAsUnmut <- length(unique(TCGA_X_As_HPV_Neg %>%
                                      filter(Patient_ID %in% pl_as) %>% .$Tumor_Sample_Barcode)) - FisherMat$TAsMut

FisherMat$InhouseMut <- sapply(g, function(x) is.integer0(IhCount[IhCount$Hugo_Symbol == x,2]))
FisherMat$InhouseUnmut <- length(unique(Inhouse_tumor_X %>%
                                          filter(Tumor_Sample_Barcode %in% pl_ih) %>% .$Tumor_Sample_Barcode)) - FisherMat$InhouseMut

FisherMat$icgcMut <- sapply(g, function(x) is.integer0(icgcCount[icgcCount$Hugo_Symbol == x,2]))
FisherMat$icgcUnmut <- n_ICGC_HPV_Neg - FisherMat$icgcMut

#combine two
FisherMat$TwoAsMut <- FisherMat$TAsMut + FisherMat$InhouseMut
FisherMat$TwoAsUnmut <- FisherMat$TAsUnmut + FisherMat$InhouseUnmut

#combine three
FisherMat$ThreeAsMut <- FisherMat$TAsMut + FisherMat$InhouseMut + FisherMat$icgcMut
FisherMat$ThreeAsUnmut <- FisherMat$TAsUnmut + FisherMat$InhouseUnmut + FisherMat$icgcUnmut

#calculate p_value
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
FisherMat$Total <- FisherMat$CauMut + FisherMat$TAsMut + FisherMat$InhouseMut + FisherMat$icgcMut
FisherMat$Gene <- rownames(FisherMat)

#write_clip(FisherMat %>% filter(p1_2As < 0.05 | p1_2As == 0.05))
#write_clip(FisherMat %>% filter(p1_3As < 0.05 | p1_3As == 0.05))
#####

#get two asian dmg
dmg3 <- FisherMat %>% filter(p1_3As < 0.05 | p1_3As == 0.05) %>% .$Gene

#GO analysis
#####
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
ego_sig <- ego[ego$pvalue < 0.05, asis =T]
barplot(ego_sig,x='Count', color = 'pvalue', showCategory = 30)
#####

#KEGG analysis
#####
gene.df <- bitr(dmg3, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)
kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") # setReadable to see the matched genes
kk2_sig <- kk2[kk2$pvalue < 0.05, asis =T]
barplot(kk2_sig,x='Count', color = 'pvalue', showCategory = 20)
#####

#trim data for plotting
trim_ca <- TCGA_X_Ca_HPV_Neg_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode != "TCGA-UF-A71A-06") %>%
  filter(Patient_ID %in% pl_ca) %>%
  filter(Hugo_Symbol %in% dmg3)
trim_as <- TCGA_X_As_HPV_Neg_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Patient_ID %in% pl_as) %>%
  filter(Hugo_Symbol %in% dmg3)
trim_ih <- Inhouse_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode %in% pl_ih) %>%
  filter(Hugo_Symbol %in% dmg3)
trim_icgc <- ICGC_HPV_Neg_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Hugo_Symbol %in% dmg3)

#Barplot all DMG dodge separate base
#####
#calculate percentage
pct <- data.frame(Gene = dmg3)
pct$`TCGA Caucasian` <- sapply(dmg3, function(x,y) {
  y = trim_ca
  a = 0 - length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_ca)
  return(a)
})
pct$`Inhouse` <- sapply(dmg3, function(x,y) {
  y = trim_ih
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_ih)
  return(a)
})
pct$`TCGA Asian` <- sapply(dmg3, function(x,y) {
  y = trim_as
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_as)
  return(a)
})
pct$`ICGC` <- sapply(dmg3, function(x,y) {
  y = trim_icgc
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / n_ICGC_HPV_Neg
  return(a)
})
#label protein size
prot <- cds %>%
  filter(feature == "CDS") %>%
  select(c(symbol, product_length)) %>%
  filter(symbol %in% dmg3) %>%
  group_by(symbol) %>%
  summarise(cds = max(product_length))
prot <- rbind(prot, c("AC138647.1", 188))
prot <- rbind(prot, c("GVINP1", 2422))
pct$Prot <- vlookup(pct$Gene, prot, "cds")
sum(is.na(pct$Prot))
pct$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
#rank gene to pvalue
row_order <- FisherMat %>% 
  filter(p1_3As < 0.05 | p1_3As == 0.05)
row_order$Prot <- vlookup(pct$Gene, prot, "cds")
row_order$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
row_order <- row_order %>%
  arrange(p1_3As) %>%
  .$label
pct <- pct %>%
  select(-c(Gene, Prot)) %>%
  gather(Cohort, Percentage, -label) %>%
  mutate(label = factor(label, levels = rev(row_order), ordered = TRUE))
max(abs(pct$Percentage))

#barplot
tiff(file = "./barplot/20200229_3as_comparison_hpv_neg_oral_2.tiff", height=25, width=6, units="in", res=100, compression="lzw")
ggplot() +
  geom_bar(data=pct[pct$Cohort=="TCGA Caucasian",], width=.75,
           stat="identity", position="stack",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_bar(data=pct[pct$Cohort!="TCGA Caucasian",], width=.75,
           stat="identity", position="dodge",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_hline(yintercept=0, lwd=0.5) +
  annotate("text", x = 2, y=0.35, label = "Asian", size=5) +
  annotate("text", x = 2, y=-0.35, label = "Caucasian", size=5) +
  
  scale_x_discrete(limits = levels(pct$Gene),expand=c(0,0)) +
  scale_fill_manual(values = c("dodgerblue3", "darkgoldenrod2", "gray50", "coral3")) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(-0.5, 0.5),
                     breaks=c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels=c("50","25","0","25","50")) +
  coord_flip() +
  xlab("Differentially mutated genes") +
  ylab("Percentage of tumors affected") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.background = element_blank()) +
  theme(axis.title.x = element_text(vjust=0,face="bold", size=12),
        axis.text.x = element_text(vjust=0, size=12)) +
  theme(axis.title.y = element_text(angle=90, vjust=1, face="bold", size=14),
        axis.text.y = element_text(size=12, face = "italic")) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.margin = unit(c(0.1,2.5,0.3,0.2),"lines"))
dev.off()
#####

#Barplot DMG > 1 dodge separate base
#####
dmg3_subset <- FisherMat %>% 
  filter(p1_3As < 0.05 | p1_3As == 0.05) %>% 
  filter(CauMut > 1 | CauMut == 1) %>%
  filter(TAsMut > 1 | TAsMut == 1) %>%
  filter(InhouseMut > 1 | InhouseMut == 1) %>%
  filter(icgcMut > 1 | icgcMut == 1) %>%
  .$Gene

pct <- data.frame(Gene = dmg3_subset)

pct$`TCGA Caucasian` <- sapply(dmg3_subset, function(x,y) {
  y = trim_ca
  a = 0 - length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_ca)
  return(a)
})

pct$`Inhouse` <- sapply(dmg3_subset, function(x,y) {
  y = trim_ih
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_ih)
  return(a)
})

pct$`TCGA Asian` <- sapply(dmg3_subset, function(x,y) {
  y = trim_as
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_as)
  return(a)
})

pct$`ICGC` <- sapply(dmg3_subset, function(x,y) {
  y = trim_icgc
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / n_ICGC_HPV_Neg
  return(a)
})
prot <- cds %>%
  filter(feature == "CDS") %>%
  select(c(symbol, product_length)) %>%
  filter(symbol %in% dmg3) %>%
  group_by(symbol) %>%
  summarise(cds = max(product_length))
prot <- rbind(prot, c("AC138647.1", 188))
prot <- rbind(prot, c("GVINP1", 2422))
pct$Prot <- vlookup(pct$Gene, prot, "cds")
sum(is.na(pct$Prot))
pct$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
#rank gene to pvalue
row_order <- FisherMat %>% 
  filter(p1_3As < 0.05 | p1_3As == 0.05) %>%
  filter(CauMut > 1 | CauMut == 1) %>%
  filter(TAsMut > 1 | TAsMut == 1) %>%
  filter(InhouseMut > 1 | InhouseMut == 1) %>%
  filter(icgcMut > 1 | icgcMut == 1)
row_order$Prot <- vlookup(pct$Gene, prot, "cds")
row_order$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
row_order <- row_order %>%
  arrange(p1_3As) %>%
  .$label
pct <- pct %>%
  select(-c(Gene, Prot)) %>%
  gather(Cohort, Percentage, -label) %>%
  mutate(label = factor(label, levels = rev(row_order), ordered = TRUE))
max(abs(pct$Percentage))

#barplot
tiff(file = "./barplot/20200229_3as_comparison_hpv_neg_oral_over1mut_2.tiff", height=5, width=6, units="in", res=100, compression="lzw")
ggplot() +
  geom_bar(data=pct[pct$Cohort=="TCGA Caucasian",], width=.75,
           stat="identity", position="stack",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_bar(data=pct[pct$Cohort!="TCGA Caucasian",], width=.75,
           stat="identity", position="dodge",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_hline(yintercept=0, lwd=0.5) +
  annotate("text", x = 1, y=0.35, label = "Asian", size=5) +
  annotate("text", x = 1, y=-0.35, label = "Caucasian", size=5) +
  
  scale_x_discrete(limits = levels(pct$label),expand=c(0,0)) +
  scale_fill_manual(values = c("dodgerblue3", "darkgoldenrod2", "gray50", "coral3")) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(-0.5, 0.5),
                     breaks=c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels=c("50","25","0","25","50")) +
  coord_flip() +
  xlab("Differentially mutated genes") +
  ylab("Percentage of tumors affected") +
  theme_bw() +
  theme(legend.position = c(0.95, 0.65)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.background = element_blank()) +
  theme(axis.title.x = element_text(vjust=0,face="bold", size=12),
        axis.text.x = element_text(vjust=0, size=12)) +
  theme(axis.title.y = element_text(angle=90, vjust=1, face="bold", size=14),
        axis.text.y = element_text(size=12, face = "italic")) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.margin = unit(c(0.1,2.5,0.3,0.2),"lines"))
dev.off()
#####

#Barplot all DMG stack, unite base
#####
#calculate percentage
pct <- data.frame(Gene = dmg3)

pct$`TCGA Caucasian` <- sapply(dmg3, function(x,y) {
  y = trim_ca
  a = 0 - length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_ca)
  return(a)
})

pct$`Inhouse` <- sapply(dmg3, function(x,y) {
  y = trim_ih
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / (length(pl_as) + length(pl_ih) + n_ICGC_HPV_Neg)
  return(a)
})

pct$`TCGA Asian` <- sapply(dmg3, function(x,y) {
  y = trim_as
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / (length(pl_as) + length(pl_ih) + n_ICGC_HPV_Neg)
  return(a)
})

pct$`ICGC` <- sapply(dmg3, function(x,y) {
  y = trim_icgc
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / (length(pl_as) + length(pl_ih) + n_ICGC_HPV_Neg)
  return(a)
})
#label protein size
prot <- cds %>%
  filter(feature == "CDS") %>%
  select(c(symbol, product_length)) %>%
  filter(symbol %in% dmg3) %>%
  group_by(symbol) %>%
  summarise(cds = max(product_length))
prot <- rbind(prot, c("AC138647.1", 188))
prot <- rbind(prot, c("GVINP1", 2422))
pct$Prot <- vlookup(pct$Gene, prot, "cds")
sum(is.na(pct$Prot))
pct$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
#rank gene to pvalue
row_order <- FisherMat %>% 
  filter(p1_3As < 0.05 | p1_3As == 0.05)
row_order$Prot <- vlookup(pct$Gene, prot, "cds")
row_order$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
row_order <- row_order %>%
  arrange(p1_3As) %>%
  .$label
pct <- pct %>%
  select(-c(Gene, Prot)) %>%
  gather(Cohort, Percentage, -label) %>%
  mutate(label = factor(label, levels = rev(row_order), ordered = TRUE))
max(abs(pct$Percentage))
#barplot
tiff(file = "./barplot/20200229_3as_comparison_hpv_neg_oral_stack_2.tiff", height=25, width=6, units="in", res=100, compression="lzw")
ggplot() +
  geom_bar(data=pct[pct$Cohort=="TCGA Caucasian",], width=.75,
           stat="identity", position="stack",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_bar(data=pct[pct$Cohort!="TCGA Caucasian",], width=.75,
           stat="identity", position="stack",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_hline(yintercept=0, lwd=0.5) +
  annotate("text", x = 5, y=0.35, label = "Asian", size=5) +
  annotate("text", x = 5, y=-0.35, label = "Caucasian", size=5) +
  
  scale_x_discrete(limits = levels(pct$label),expand=c(0,0)) +
  scale_fill_manual(values = c("dodgerblue3", "darkgoldenrod2", "gray50", "coral3")) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(-0.5, 0.5),
                     breaks=c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels=c("50","25","0","25","50")) +
  coord_flip() +
  xlab("Differentially mutated genes") +
  ylab("Percentage of tumors affected") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.background = element_blank()) +
  theme(axis.title.x = element_text(vjust=0,face="bold", size=12),
        axis.text.x = element_text(vjust=0, size=12)) +
  theme(axis.title.y = element_text(angle=90, vjust=1, face="bold", size=14),
        axis.text.y = element_text(size=12, face = "italic")) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.margin = unit(c(0.1,2.5,0.3,0.2),"lines"))
dev.off()
#####

#Barplot DMG > 1 stack, unite base
#####
dmg3_subset <- FisherMat %>% 
  filter(p1_3As < 0.05 | p1_3As == 0.05) %>% 
  filter(CauMut > 1 | CauMut == 1) %>%
  filter(TAsMut > 1 | TAsMut == 1) %>%
  filter(InhouseMut > 1 | InhouseMut == 1) %>%
  filter(icgcMut > 1 | icgcMut == 1) %>%
  .$Gene

pct <- data.frame(Gene = dmg3_subset)

pct$`TCGA Caucasian` <- sapply(dmg3_subset, function(x,y) {
  y = trim_ca
  a = 0 - length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / length(pl_ca)
  return(a)
})

pct$`Inhouse` <- sapply(dmg3_subset, function(x,y) {
  y = trim_ih
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / (length(pl_as) + length(pl_ih) + n_ICGC_HPV_Neg)
  return(a)
})

pct$`TCGA Asian` <- sapply(dmg3_subset, function(x,y) {
  y = trim_as
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / (length(pl_as) + length(pl_ih) + n_ICGC_HPV_Neg)
  return(a)
})

pct$`ICGC` <- sapply(dmg3_subset, function(x,y) {
  y = trim_icgc
  a = length(unique(filter(y, Hugo_Symbol == x)$Tumor_Sample_Barcode)) / (length(pl_as) + length(pl_ih) + n_ICGC_HPV_Neg)
  return(a)
})
#label protein size
prot <- cds %>%
  filter(feature == "CDS") %>%
  select(c(symbol, product_length)) %>%
  filter(symbol %in% dmg3) %>%
  group_by(symbol) %>%
  summarise(cds = max(product_length))
prot <- rbind(prot, c("AC138647.1", 188))
prot <- rbind(prot, c("GVINP1", 2422))
pct$Prot <- vlookup(pct$Gene, prot, "cds")
sum(is.na(pct$Prot))
pct$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
#rank gene to pvalue
row_order <- FisherMat %>% 
  filter(p1_3As < 0.05 | p1_3As == 0.05) %>%
  filter(CauMut > 1 | CauMut == 1) %>%
  filter(TAsMut > 1 | TAsMut == 1) %>%
  filter(InhouseMut > 1 | InhouseMut == 1) %>%
  filter(icgcMut > 1 | icgcMut == 1)
row_order$Prot <- vlookup(pct$Gene, prot, "cds")
row_order$label <- paste(pct$Gene, " (", pct$Prot, "aa)", sep = "")
row_order <- row_order %>%
  arrange(p1_3As) %>%
  .$label
pct <- pct %>%
  select(-c(Gene, Prot)) %>%
  gather(Cohort, Percentage, -label) %>%
  mutate(label = factor(label, levels = rev(row_order), ordered = TRUE))
max(abs(pct$Percentage))

#barplot
tiff(file = "./barplot/20200229_3as_comparison_hpv_neg_oral_over1mut_stack_2.tiff", height=5, width=5, units="in", res=100, compression="lzw")
ggplot() +
  geom_bar(data=pct[pct$Cohort=="TCGA Caucasian",], width=.75,
           stat="identity", position="stack",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_bar(data=pct[pct$Cohort!="TCGA Caucasian",], width=.75,
           stat="identity", position="stack",
           aes(x=label, y=Percentage, fill=Cohort), color = "black") +
  geom_hline(yintercept=0, lwd=0.5) +
  annotate("text", x = 1, y=0.3, label = "Asian", size=5) +
  annotate("text", x = 1, y=-0.3, label = "Caucasian", size=5) +
  
  scale_x_discrete(limits = levels(pct$label),expand=c(0,0)) +
  scale_fill_manual(values = c("dodgerblue3", "darkgoldenrod2", "gray50", "coral3")) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(-0.5, 0.5),
                     breaks=c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels=c("50","25","0","25","50")) +
  coord_flip() +
  xlab("Differentially mutated genes") +
  ylab("Percentage of tumors affected") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.7)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.background = element_blank()) +
  theme(axis.title.x = element_text(vjust=0,face="bold", size=12),
        axis.text.x = element_text(vjust=0, size=12)) +
  theme(axis.title.y = element_text(angle=90, vjust=1, face="bold", size=14),
        axis.text.y = element_text(size=12, face = "italic")) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.margin = unit(c(0.1,2.5,0.3,0.2),"lines"))
dev.off()
#####