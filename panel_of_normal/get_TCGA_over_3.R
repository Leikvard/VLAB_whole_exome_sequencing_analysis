library(tidyverse)

#TCGA_all_mut is obtained from another R script, in tcga_33_cancer_type_maf folder
TCGA_all_mut <- TCGA_all_mut %>% unite("Lift_Pos", c("Chromosome", "Start_Position"), sep = ":", remove = FALSE)
TCGA_all_mut <- TCGA_all_mut %>% unite("Lift_Pos2", c("Lift_Pos", "End_Position"), sep = "-", remove = FALSE)
TCGA_all_mut <- TCGA_all_mut %>% unite("Check_count", c("Lift_Pos2", "Reference_Allele","Tumor_Seq_Allele2"), sep = "-", remove = FALSE)
TCGA3 <- TCGA_all_mut %>% count(Check_count)
TCGA3 <- TCGA3 %>% filter(n > 3 | n == 3)
TCGA3 <- unique(filter(TCGA_all_mut, Check_count %in% TCGA3$Check_count)$Lift_Pos2)
write.csv(TCGA3, "TCGA_over_3_hg38.csv")

x <- TCGA_all_mut %>%
  select(c(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode)) %>%
  unite("N1", c("Chromosome", "Start_Position"), sep = ":", remove = "TRUE", na.rm = FALSE) %>%
  unite("N2", c("N1", "End_Position"), sep = "-", remove = "TRUE", na.rm = FALSE) %>%
  unite("N3", c("N2", "Reference_Allele"), sep = "", remove = "FALSE", na.rm = FALSE) %>%
  unite("Genomic_Change", c("N3", "Tumor_Seq_Allele2"), sep = ">", remove = "FALSE", na.rm = FALSE)

y <- x %>%
  group_by(Genomic_Change) %>% 
  summarise(n = length(unique(Tumor_Sample_Barcode))) %>%
  filter(n > 3 | n == 3)

TCGA_hg38_over3 <- x %>%
  filter(Genomic_Change %in% y$Genomic_Change)

TCGA_hg38_over3 <- unique(TCGA_hg38_over3$N2)

write.csv(TCGA_hg38_over3, "TCGA_over_3_hg38.csv")

