library(readxl)
library(tidyverse)
TypeDef <- read_excel("NONSYN SYN RNA mutations DNA link consensus_VL_19Dec2019.xlsx", sheet = 2)
NS <- as.vector(TypeDef$`Non-syn`)
#Inhouse_cell_X$Hugo_Symbol <- str_remove(Inhouse_cell_X$Hugo_Symbol, "`")

gene_list <- read_excel("Immune_gene_list.xlsx")
gene_list <- gene_list$SYMBOL


> count_tumor_total <- Inhouse_tumor_X %>% filter(Variant_Classification %in% NS) %>% count(ID)
> count_tumor_immune <- Inhouse_tumor_X %>% filter(Variant_Classification %in% NS) %>% filter(Hugo_Symbol %in% gene_list) %>% count(ID)
> df_inhouse_tumor <- data.frame(ID = count_tumor_total$ID, Immune = count_tumor_immune$n, Total = count_tumor_total$n)
> df_inhouse_tumor$Percentage <- df_inhouse_tumor$Immune / df$Total


> count_cell_total <- Inhouse_cell_X %>% filter(Variant_Classification %in% NS) %>% count(ID)
> count_cell_immune <- Inhouse_cell_X %>% filter(Variant_Classification %in% NS) %>% filter(Hugo_Symbol %in% gene_list) %>% count(ID)
> df_inhouse_cell <- data.frame(ID = count_cell_total$ID, Immune = count_cell_immune$n, Total = count_cell_total$n)
> df_inhouse_cell$Percentage <- df_inhouse_cell$Immune / df_inhouse_cell$Total


count_tcga_Asian_total <- TCGA_X_As %>% filter(Variant_Classification %in% NS) %>% count(Patient_ID)
> count_tcga_Asian_immune <- TCGA_X_As %>% filter(Variant_Classification %in% NS) %>% filter(Hugo_Symbol %in% gene_list) %>% count(Patient_ID)
> df_asian <- data.frame(ID = count_tcga_Asian_total$Patient_ID, Immune = count_tcga_Asian_immune$n, Total = count_tcga_Asian_total$n)
> df_asian$Percentage <- df_asian$Immune / df_asian$Total


> count_tcga_caucasian_total <- TCGA_X_Ca %>% filter(Variant_Classification %in% NS) %>% count(Patient_ID)
> count_tcga_caucasian_immune <- TCGA_X_Ca %>% filter(Variant_Classification %in% NS) %>% filter(Hugo_Symbol %in% gene_list) %>% count(Patient_ID)
> which(count_tcga_caucasian_total$Patient_ID %in% count_tcga_caucasian_immune$Patient_ID == FALSE)
[1]  19 126 279 331
> df_caucasian <- data.frame(ID = count_tcga_caucasian_total$Patient_ID[-c(19, 126, 279, 331)], Immune = count_tcga_caucasian_immune$n, Total = count_tcga_caucasian_total$n[-c(19, 126, 279, 331)])
> df_caucasian$Percentage <- df_caucasian$Immune / df_caucasian$Total
> df_caucasian$ID <- as.character(df_caucasian$ID)
> df_caucasian[433,] <- c(count_tcga_caucasian_total$Patient_ID[19],0,count_tcga_caucasian_total$n[19],0)
> df_caucasian[434,] <- c(count_tcga_caucasian_total$Patient_ID[126],0,count_tcga_caucasian_total$n[126],0)
> df_caucasian[435,] <- c(count_tcga_caucasian_total$Patient_ID[279],0,count_tcga_caucasian_total$n[279],0)
> df_caucasian[436,] <- c(count_tcga_caucasian_total$Patient_ID[331],0,count_tcga_caucasian_total$n[331],0)



> count_tcga_caucasian_hpv_pos_total <- TCGA_X_Ca_HPV_Pos %>% filter(Variant_Classification %in% NS) %>% count(Patient_ID)
> count_tcga_caucasian_hpv_pos_immune <- TCGA_X_Ca_HPV_Pos %>% filter(Variant_Classification %in% NS) %>% filter(Hugo_Symbol %in% gene_list) %>% count(Patient_ID)
> df_caucasian$Percentage <- df_caucasian$Immune / df_caucasian$Total
> df_caucasian$ID <- as.character(df_caucasian$ID)
> which(count_tcga_caucasian_total$Patient_ID %in% count_tcga_caucasian_immune$Patient_ID == FALSE)
> which(count_tcga_caucasian_hpv_pos_total$Patient_ID %in% count_tcga_caucasian_hpv_pos_immune$Patient_ID == FALSE)
> df_caucasian_hpv_pos[71,] <- c(count_tcga_caucasian_hpv_pos_total$Patient_ID[5],0,count_tcga_caucasian_hpv_pos_total$n[5],0)
> df_caucasian_hpv_pos[72,] <- c(count_tcga_caucasian_hpv_pos_total$Patient_ID[49],0,count_tcga_caucasian_hpv_pos_total$n[49],0)

count_tcga_caucasian_hpv_neg_total <- TCGA_X_Ca_HPV_Neg %>% filter(Variant_Classification %in% NS) %>% count(Patient_ID)
count_tcga_caucasian_hpv_neg_immune <- TCGA_X_Ca_HPV_Neg %>% filter(Variant_Classification %in% NS) %>% filter(Hugo_Symbol %in% gene_list) %>% count(Patient_ID)
which(count_tcga_caucasian_hpv_neg_total$Patient_ID %in% count_tcga_caucasian_hpv_neg_immune$Patient_ID == FALSE)
df_caucasian_hpv_neg <- data.frame(ID = count_tcga_caucasian_hpv_neg_total$Patient_ID[-c(97,222)], Immune = count_tcga_caucasian_hpv_neg_immune$n, Total = count_tcga_caucasian_hpv_neg_total$n[-c(97,222)])
df_caucasian_hpv_neg$Percentage <- df_caucasian_hpv_neg$Immune / df_caucasian_hpv_neg$Total
df_caucasian_hpv_neg$ID <- as.character(df_caucasian_hpv_neg$ID)
df_caucasian_hpv_neg[348,] <- c(count_tcga_caucasian_hpv_neg_total$Patient_ID[97],0,count_tcga_caucasian_hpv_neg_total$n[97],0)
df_caucasian_hpv_neg[349,] <- c(count_tcga_caucasian_hpv_neg_total$Patient_ID[222],0,count_tcga_caucasian_hpv_neg_total$n[222],0)

