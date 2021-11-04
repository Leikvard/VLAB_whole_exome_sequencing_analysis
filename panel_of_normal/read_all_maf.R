#In terminal
#move to top tree directory
#Run in terminal: gzip -d $(find ./ -type f -name '*.gz')
#Run in terminal: mv **/*.maf /path/to/foler/want/to/move
#To gzip back find . \( -name '*.maf' -o -name '*.js' \) -exec gzip --verbose --keep --best --force {} \;

library(tidyverse)
require(data.table)
file_list <- list.files()
file_list <- file_list[-1]
file_list <- file_list[-1]
TCGA_all_mut <- data.frame()

for (i in 1:length(file_list)){
  temp_data <- fread(file_list[i], stringsAsFactors = F) #read in files using the fread function from the data.table package
  TCGA_all_mut <- rbindlist(list(TCGA_all_mut, temp_data), use.names = T) #for each iteration, bind the new data to the building TCGA_all_mut
}
#TCGA_all_mut <- TCGA_all_mut %>% select(c(1,5,6,7,9,16,36,37))
