#In terminal
#move to top tree directory
#Run in terminal: gzip -d $(find ./ -type f -name '*.gz')
#Run in terminal: mv **/*.maf /path/to/foler/want/to/move
#To gzip back find . \( -name '*.maf' -o -name '*.js' \) -exec gzip --verbose --keep --best --force {} \;

library(tidyverse)
require(data.table)
file_list <- list.files()
file_list <- file_list[-4]
Inhouse_all_blood <- data.frame()

for (i in 1:length(file_list)){
  temp_data <- fread(file_list[i], stringsAsFactors = F) #read in files using the fread function from the data.table package
  Inhouse_all_blood <- rbindlist(list(Inhouse_all_blood, temp_data), use.names = F) #for each iteration, bind the new data to the building Inhouse_all_blood
}
#Inhouse_all_blood <- Inhouse_all_blood %>% select(c(1,5,6,7,9,16,36,37))

Export <- Inhouse_all_blood %>%
  unite("Genomic_position", c(chr_location, POS), sep = "-", remove = FALSE) %>%
  select(c(Genomic_position, chr_location))

write.csv(Export, file = "inhouse_blood_all_mut.csv")




