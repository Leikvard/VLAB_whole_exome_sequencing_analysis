library(tidyverse)
#Read COSMIC variants
require(data.table)
COSMIC <- data.table::fread("CosmicMutantExport.tsv")
xyz <- COSMIC %>% select(c(GENOMIC_MUTATION_ID, ID_tumour)) %>%
  group_by(GENOMIC_MUTATION_ID) %>% 
  summarise(n = length(unique(ID_tumour))) %>%
  filter(n > 10 | n == 10) %>%
  filter(GENOMIC_MUTATION_ID != "")

COSMIC_pos <- COSMIC %>% 
  filter(GENOMIC_MUTATION_ID %in% xyz$GENOMIC_MUTATION_ID) %>%
  select(`Mutation genome position`) %>%
  separate(`Mutation genome position`, into = c("P1", "P2"), remove = TRUE, sep = "-")

COSMIC_pos <- unique(COSMIC_pos$P1)

write.csv(COSMIC_pos, "COSMIC_genomic_position_over_10.csv")
