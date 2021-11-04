library(tidyverse)
library(expss)

QM36L <- read.csv("./Environment/QM36_line_filtered_20200221.csv")
gene.df <- bitr("KAT14", fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Cell_copy <- Inhouse_cell_X
x <- as.data.frame(matrix(, nrow = nrow(QM36L), ncol = ncol(Cell_copy)))
colnames(x) <- colnames(Cell_copy)
x$Hugo_Symbol <- QM36L$Gene_Name
x$Tumor_Sample_Barcode <- "QM36"
x$Protein_Change <- QM36L$HGVS.p
x$Chromosome <- QM36L$CHROM
x$Start_position <- QM36L$POS
x$Reference_Allele <- QM36L$REF
x$Tumor_Seq_Allele2 <- QM36L$ALT
x$dbSNP_RS <- QM36L$ID
x$Variant_Type <- QM36L$Variant.Type
x$AF <- QM36L$W095_AD.ratio.Alt.Ref.Alt.
x$Variant_Classification <- paste(QM36L$Variant.Type, QM36L$Effect)
x$Variant_Classification <- gsub("SNP missense_variant", "Missense_Mutation", x$Variant_Classification)
x$Variant_Classification <- gsub("SNP stop_gained", "Nonsense_Mutation", x$Variant_Classification)
x$Variant_Classification <- gsub("DEL,INS inframe_insertion", "In_Frame_Ins", x$Variant_Classification)
x$Variant_Classification <- gsub("INS frameshift_variant", "Frame_Shift_Ins", x$Variant_Classification)
x$Variant_Classification <- gsub("DEL frameshift_variant", "Frame_Shift_Del", x$Variant_Classification)
x$Variant_Classification <- gsub("INS inframe_insertion", "In_Frame_Ins", x$Variant_Classification)
x$Variant_Classification <- gsub("SNP,INTERVAL missense_variant", "Missense_Mutation", x$Variant_Classification)
x$Variant_Classification <- gsub("SNP splice_region_variant&synonymous_variant", "Silent", x$Variant_Classification)
x$Variant_Classification <- gsub("&splice_region_variant", "", x$Variant_Classification)
x$Variant_Classification <- gsub("INS disruptive_inframe_insertion", "In_Frame_Ins", x$Variant_Classification)
x$Variant_Classification <- gsub("INS,In_Frame_Ins", "In_Frame_Ins", x$Variant_Classification)
x$Variant_Classification <- gsub("DEL disruptive_inframe_deletion", "In_Frame_Del", x$Variant_Classification)
x$Variant_Classification <- gsub("DEL,INS disruptive_inframe_deletion", "In_Frame_Del", x$Variant_Classification)

Cell_copy <- rbind(Cell_copy, x)
Cell_copy2 <- Cell_copy
Cell_copy <- as.data.frame(Cell_copy)
  write_clip(unique(Cell_copy$Hugo_Symbol))


#convert inhouse
x <- read.csv("./Environment/cell_hugo.csv", header = FALSE)
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
x <- x[-which(x$Input == "MUM1" & x$Hugo != "PWWP3A"),] #remove MUM1 that is translated to IRF4
x <- x[-which(x$Input == "LOR"),]
x <- x[-which(x$Input == "DEC1" & x$Hugo != "DELEC1"),]
x <- x[-which(x$Input == "CSRP2BP" & x$Hugo != "KAT14"),]
x$Hugo[duplicated(x$Hugo)]
#CHDC2 and CXorf22 are both converted to CFAP47, they are all on X, position close to each other, probably same gene
y$Hugo[duplicated(y$Hugo)]
dic <- rbind(x,y)
dic$Input[duplicated(dic$Input)]
dic$Hugo[duplicated(dic$Hugo)]
#all checked, inputs with same hugo are actually same genes
Cell_copy$Hugo_Symbol <- sapply(Cell_copy2$Hugo_Symbol, function(x) {
  if (x %in% dic$Input) {
    return(vlookup(x, dic, "Hugo"))
  }
  else {
    return(x)
  }
})
Cell_copy <- as.data.frame(Cell_copy)
colnames(Cell_copy)[1] <- "Cell_line"

gene_list <- read.csv("./gene_list/Genelist.csv")
gene_list <- gene_list %>% filter(Pathway == "Drug") %>% .$Gene
gene_list <- as.character(gene_list)
ID = "QM13"
x <- Cell_copy %>% 
  filter(Hugo_Symbol %in% cancer_gene_census_combined) %>%
  filter(Cell_line == ID) %>%
  filter(Variant_Classification %in% NS) %>%
  select(c(Cell_line, Hugo_Symbol, AF, Protein_Change)) %>%
  arrange(desc(AF))
write_clip(x)

x <- Cell_copy %>% 
  filter(Variant_Classification %in% NS) %>%
  filter(Cell_line != "QM25-P04") %>%
  filter(Cell_line != "QM25-P82") %>%
  filter(Cell_line != "T48-P08") %>%
  filter(Cell_line != "T48-P95") %>%
  filter(Cell_line != "T73-PD") %>%
  filter(Hugo_Symbol %in% gene_list) %>%
  group_by(Hugo_Symbol) %>%
  summarise(n = length(unique(Cell_line))) %>%
  arrange(desc(n))


write_clip(x)
write_clip(NS)
ID <- c("QM13", "QM16", "QM17", "QM25", 'QM28',
        'QM29','QM43',
        'T48','T63','T76','T77')
df <- as.data.frame(matrix(, nrow = 0, ncol = 4))
colnames(df) <- c("Cell_line", 'Hugo_Symbol', 'AF', "Protein_Change")
for (i in c(1:11)){
  id = ID[i]
  x <- Cell_copy %>% 
    filter(Hugo_Symbol %in% gene_list) %>%
    filter(Cell_line == id) %>%
    filter(Variant_Classification %in% NS) %>%
    select(c(Cell_line, Hugo_Symbol, AF, Protein_Change)) %>%
    arrange(desc(AF))
  df <- rbind(df,x)
}
write_clip(df)

ID = "QM36"
x <- Cell_copy %>% 
  filter(Hugo_Symbol %in% gene_list) %>%
  filter(Cell_line == ID) %>%
  select(c(Cell_line, Hugo_Symbol, AF, Protein_Change, Start_position, )) %>%
  arrange(desc(AF))
write_clip(x)
