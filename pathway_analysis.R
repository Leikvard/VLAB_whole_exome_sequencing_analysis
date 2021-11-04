library(tidyverse)
library(readxl)
library(ggpubr)

#Read in pathway definitions
Pathway_List <- read.csv(file = "./gene_list/Genelist.csv")
Pathway_Name <- as.vector(unique(Pathway_List$Pathway))

#hpv-negative oscc tumors
#trim data sets
#####
#Caucasian
pl_ca <- TCGA_Pat_Ca_HPV_Pos %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_ca <- as.character(pl_ca)
trim_ca <- TCGA_X_Ca_HPV_Pos_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Patient_ID %in% pl_ca)
#icgc hpv-pos ns
trim_ic <- ICGC_HPV_Pos_copy %>% 
  filter(Variant_Classification %in% NS)
#####

# Construct result dataframe
#####
df <- matrix(0, nrow = 7, ncol = 4)
colnames(df) <- c("TCGA_Caucasian_mutated", "TCGA_Caucasian_unmutated", 
                  "ICGC_mutated", "ICGC_unmutated")
rownames(df) <- Pathway_Name[1:7]

#Read in result dataframe
for (i in c(1:7)){
  p <- Pathway_Name[i]
  g <- as.vector(filter(Pathway_List, Pathway == p)$Gene)
  tc <- filter(trim_ca, Hugo_Symbol %in% g)
  ti <- filter(trim_ic, Hugo_Symbol %in% g)
  
  df[p,1] <- length(unique(tc$Tumor_Sample_Barcode))
  df[p,2] <- length(unique(trim_ca$Tumor_Sample_Barcode)) - df[p,1]
  df[p,3] <- length(unique(ti$Tumor_Sample_Barcode))
  df[p,4] <- length(unique(trim_ic$Tumor_Sample_Barcode)) - df[p,3]
}
df <- as.data.frame(df)
#Calculate P value and Odds ratio
df$p1_2As <- 
  apply(df,1,function(x)
    if (fisher.test(matrix(x[c(1:4)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "less")$p.value)
    }
    else {
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "greater")$p.value)
    }
  )

df$odds_2As <- 
  apply(df,1,function(x)
    if (fisher.test(matrix(x[c(1:4)],nrow=2))$estimate < 1){
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "less")$estimate)
    }
    else {
      return(fisher.test(matrix(x[c(1:4)],nrow=2), alternative = "greater")$estimate)
    }
  )


write.csv(df, "./result_table/pathway_analysis_positive.csv")
#####


#Plotting Bar graph, dodge, separate base
#####
#pathway df
Bdf <- df %>% 
  select(c("TCGA_Caucasian_mutated", "ICGC_mutated")) 
Bdf$TCGA_Caucasian_mutated <- 0 - Bdf$TCGA_Caucasian_mutated / length(unique(trim_ca$Tumor_Sample_Barcode))
Bdf$ICGC_mutated <- Bdf$ICGC_mutated / n_ICGC_HPV_Pos
colnames(Bdf) <- c("TCGA Caucasian tumor", "ICGC tumor")
Bdf$Pathway <- rownames(Bdf)
Bdf <- Bdf %>% gather(Cohort, Percentage, -Pathway)
Order <- Bdf %>% 
  filter(Cohort=="TCGA Caucasian tumor") %>% 
  arrange(Percentage) %>% .$Pathway
Bdf <- Bdf %>%
  mutate(Pathway=factor(Pathway, levels=Order, ordered=TRUE))

#p-value df
Bdf1 <- df %>% select(p1_2As)
Bdf1$p1 <- -log10(Bdf1$p1)
Bdf1$Pathway <- rownames(Bdf1)
Bdf1 <- Bdf1 %>% gather(Cohort, p_value, -Pathway)
Bdf1$Cohort <- as.character(Bdf1$Cohort)
Bdf1 <- Bdf1 %>%
  mutate(Pathway=factor(Pathway, levels=Order, ordered=TRUE))

#plot
tiff(file = "./barplot/pathway_analysis_hpv_pos_oscc.tiff", height=3, width=8, units="in", res=100, compression="lzw")
p1 <- ggplot() +
  geom_bar(data=Bdf, width=.75,
           stat="identity", position="stack",
           aes(x=Pathway, y=Percentage, fill=Cohort)) +
  geom_hline(yintercept=0, lwd=1) +
  annotate("text", x = 1.5, y=0.5, label = "Asian", size=6, fontface = "bold") +
  annotate("text", x = 1.5, y=-0.5, label = "Caucasian", size=6, fontface = "bold") +
  
  scale_x_discrete(limits = rev(levels(Bdf$Pathway)),expand=c(0,0)) +
  scale_fill_manual(values = c("dodgerblue3", "coral3"),
                    drop=FALSE) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(-1, 1),
                     breaks=c(-1,-0.5,0,0.5,1),
                     labels=c("100","50","0","50","100")) +
  coord_flip() +
  xlab("Pathways") +
  ylab("Percentage of patients affected") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 8)) +
  theme(legend.background = element_blank()) +
  theme(axis.title.x = element_text(vjust=0.5,face="bold", size=10),
        axis.text.x = element_text(vjust=-1, size=10)) +
  theme(axis.title.y = element_text(angle=90, vjust=1, face="bold", size=10),
        axis.text.y = element_text(size=10)) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(size=1, color = "black")) +
  theme(plot.margin = unit(c(0.2,0.9,0.3,0.2),"lines"))

p2 <- ggplot() +
  geom_bar(data=Bdf1, width=.75,
           stat="identity", position="dodge",
           aes(x=Pathway, y=p_value, fill = p_value)) +
  geom_hline(yintercept=1.30103, lwd=0.3, lty = "dashed") +
  #annotate("text", x = 1.5, y=0.5, label = "Asian", size=8, fontface = "bold") +
  #annotate("text", x = 1.5, y=-0.5, label = "Caucasian", size=8, fontface = "bold") +
  #annotate("text", x = 3, y=0.7, label = "HPV-negative tumors", size=8, fontface = "bold") +
  
  scale_x_discrete(limits = rev(levels(Bdf1$Pathway)),expand=c(0,0)) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_y_continuous(expand = c(0,0), 
                     limits=c(0, 6)) +
  coord_flip() +
  xlab("") +
  ylab("-log(p-value)") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.55)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 8)) +
  theme(legend.background = element_blank()) +
  theme(axis.title.x = element_text(vjust=0.5,face="bold", size=10),
        axis.text.x = element_text(vjust=-1, size=10)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid=element_blank()) +
  theme(panel.border = element_rect(size=1, color = "black")) +
  theme(plot.margin = unit(c(0.2,0.9,0.3,0.2),"lines"))

figure <- ggarrange(p1, p2, align = "h",widths = c(3,1),
                    ncol = 2, nrow = 1)

annotate_figure(figure,
                top = text_grob("HPV-positive OSCC", face = "bold", size = 20, vjust = 0.3)
                #bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                #                   hjust = 1, x = 1, face = "italic", size = 10),
                #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                #right = "I'm done, thanks :-)!",
                #fig.lab = "Figure 1", fig.lab.face = "bold"
)

dev.off()
#####

