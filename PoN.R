library(tidyverse)

#Ex is all variant with AF >= 1e-5 in ExAC
#inhouse_blood is all variants occured in our inhouse patients' blood
#tcga3 is tcga variants with occurence over 3.
#cosmic10 is cosmic variants with occurrence over 10
Ex <- Ex[!duplicated(Ex$Genomic_Change),]

inhouse_blood <- read.csv("inhouse_blood_all_mut.csv")
inhouse_blood <- inhouse_blood[!duplicated(inhouse_blood$chr_location),]

tcga3 <- as.data.frame(read.table("TCGA_over_3_hg19.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
tcga3 <- tcga3 %>%
  separate(V1, into = c("N1", "N2"), remove = TRUE, sep = "-")
tcga3$N1 <- str_remove(tcga3$N1, "chr")
tcga3 <- tcga3[!duplicated(tcga3$N1),]

#tcga_all <- read.table("TCGA_all_mutation_hg19_from_json_file.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#tcga_all$V1 <- str_remove(tcga_all$V1, "chr")

cosmic10 <- read.csv("COSMIC_genomic_position_over_10.csv")

PoN <- unique(Ex$Genomic_Position)
PoN <- PoN[!PoN %in% tcga3$N1]
PoN <- PoN[!PoN %in% cosmic10$x]




#QM36T 
QM36T <- read.csv("W094_QM36_T_minus_hg19_DupRemoved_20200109.txt", sep = "\t")

QM36T$Gene_Name <- str_remove(QM36T$Gene_Name,"`")
QM36T$HGVS.p <- str_remove(QM36T$HGVS.p,"`")
QM36T$HGVS.p <- gsub("\\*", "Ter", QM36T$HGVS.p)
QM36T <- QM36T %>% unite("Query", c("Gene_Name","HGVS.p"), sep = "-", remove = FALSE)
QM36T <- QM36T[!duplicated(QM36T$Query),]

QM36T <- QM36T %>% filter(!chr_location %in% PoN)
QM36T <- QM36T %>% filter(HGVS.p != "") %>% filter(Effect != "synonymous_variant")

#length(QM36T$Query[QM36T$Query %in% TCGA_all_mut$Query])
#length(QM36T$chr_location[!QM36T$chr_location %in% inhouse_blood$chr_location])
QM36T <- QM36T %>% 
  filter(!chr_location %in% inhouse_blood$chr_location)# %>% 
  #filter(Query %in% TCGA_all_mut$Query) %>% 


#QM36L
QM36L <- read.csv("W095_QM36_P30_minus_hg19_DupRemoved_20200109.txt", sep = "\t")

QM36L$Gene_Name <- str_remove(QM36L$Gene_Name,"`")
QM36L$HGVS.p <- str_remove(QM36L$HGVS.p,"`")
QM36L$HGVS.p <- gsub("\\*", "Ter", QM36L$HGVS.p)
QM36L <- QM36L %>% unite("Query", c("Gene_Name","HGVS.p"), sep = "-", remove = FALSE)
QM36L <- QM36L[!duplicated(QM36L$Query),]

QM36L <- QM36L %>% filter(!chr_location %in% PoN)
QM36L <- QM36L %>% filter(HGVS.p != "") %>% filter(Effect != "synonymous_variant")

length(QM36L$Query[QM36L$Query %in% TCGA_all_mut$Query])
length(QM36L$chr_location[!QM36L$chr_location %in% inhouse_blood$chr_location])
nrow(QM36L[(!QM36L$chr_location %in% inhouse_blood$chr_location) & (QM36L$Query %in% TCGA_all_mut$Query),])
QM36L <- QM36L %>% 
  filter(!chr_location %in% inhouse_blood$chr_location) #%>% 
  #filter(Query %in% TCGA_all_mut$Query) %>% 


#HSC6
HSC6 <- read.csv("W027_HSC-6_P_minus_hg19_DupRemoved_20200109.txt", sep = "\t")
HSC6$Gene_Name <- str_remove(HSC6$Gene_Name,"`")
HSC6$HGVS.p <- str_remove(HSC6$HGVS.p,"`")
HSC6$HGVS.p <- gsub("\\*", "Ter", HSC6$HGVS.p)
HSC6 <- HSC6 %>% unite("Query", c("Gene_Name","HGVS.p"), sep = "-", remove = FALSE)
HSC6 <- HSC6[!duplicated(HSC6$Query),]

HSC6 <- HSC6 %>% filter(!chr_location %in% PoN)
HSC6 <- HSC6 %>% filter(HGVS.p != "") %>% filter(Effect != "synonymous_variant")

length(HSC6$Query[HSC6$Query %in% TCGA_all_mut$Query])
length(HSC6$chr_location[!HSC6$chr_location %in% inhouse_blood$chr_location])
HSC6 <- HSC6 %>% 
  filter(!chr_location %in% inhouse_blood$chr_location) #%>% 
  #filter(Query %in% TCGA_all_mut$Query) %>% 
  #.$Query


write.csv(QM36L, "QM36_line_filtered_20200221.csv")
write.csv(QM36T, "QM36_tumor_filtered_20200221.csv")
write.csv(HSC6, "HSC6_filtered_20200221.csv")


#QM29(positive control)
QM29 <- read.csv("W080_QM29_P31_minus_hg19 copy.txt", sep = "\t")
QM29 <- QM29 %>% unite(chr_location, c(CHROM, POS), sep = ":", remove = FALSE)
QM29$Gene_Name <- str_remove(QM29$Gene_Name,"`")
QM29$HGVS.p <- str_remove(QM29$HGVS.p,"`")
QM29$HGVS.p <- gsub("\\*", "Ter", QM29$HGVS.p)
QM29 <- QM29 %>% unite("Query", c("Gene_Name","HGVS.p"), sep = "-", remove = FALSE)
QM29 <- QM29[!duplicated(QM29$Query),]

QM29 <- QM29 %>% filter(!chr_location %in% PoN2)
#QM29 <- QM29 %>% filter(HGVS.p != "") %>% filter(Effect != "synonymous_variant")

#length(QM29$Query[QM29$Query %in% TCGA_all_mut$Query])
QM29 <- QM29 %>% 
  filter(!chr_location %in% inhouse_blood$chr_location)

QM29_somatic <- read_excel("All Primary cell lines vs blood_WES_DNA Link_20191211.xlsx")
QM29_somatic <- QM29_somatic %>% filter(ID == "QM29-P31") %>% 
  unite(chr_location, c(Chromosome, Start_position), sep = ":", remove = FALSE) %>%
  unite(chr_location2, c(Chromosome, End_position), sep = ":", remove = FALSE)
QM29_somatic$Hugo_Symbol <- str_remove(QM29_somatic$Hugo_Symbol, "`")
x <- QM29$chr_location
x <- x[!duplicated(x)]
y <- QM29_somatic$chr_location
z <- QM29_somatic$chr_location2
length(x[x%in%y])
length(x[x%in%z])
length(y[y%in%z])
length(y[y%in%PoN2])
QM29_unoverlap <- QM29_somatic %>% filter(chr_location %in% y[!y%in%x])


