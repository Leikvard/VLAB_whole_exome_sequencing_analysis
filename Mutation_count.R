library(tidyverse)
library(ggpubr)

#Define sample list for each cohort
pl_ca <- TCGA_Pat_Ca_HPV_Neg %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_ca <- as.character(pl_ca)
trim_ca <- TCGA_X_Ca_HPV_Neg_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode != "TCGA-UF-A71A-06") %>%
  filter(Patient_ID %in% pl_ca)

#Asian
pl_as <- TCGA_Pat_As_HPV_Neg %>% filter(MAJ_SITE == "Oral Cavity") %>% 
  .$PATIENT_ID # get patient list of caucasian hpv-negative oral cancer
pl_as <- as.character(pl_as)
trim_as <- TCGA_X_As_HPV_Neg_copy %>%
  filter(Variant_Classification %in% NS) %>%
  filter(Patient_ID %in% pl_as) 

#icgc hpv-neg ns
pl_ic <- ICGC_Pat %>% filter(HPV.status == "Neg") %>% 
  .$PATIENT_ID
trim_ic <- ICGC_HPV_Neg_copy %>% 
  filter(Variant_Classification %in% NS)

#inhouse
pl_ih <- Inhouse_Pat %>% filter(MAJ_SITE == "Oral Cavity") %>% .$PATIENT_ID
pl_ih <- as.character(pl_ih)
trim_ih <- Inhouse_copy %>% 
  filter(Variant_Classification %in% NS) %>%
  filter(Tumor_Sample_Barcode %in% pl_ih)

#Read mutation count
df <- data.frame()
temp <- data.frame(row.names = pl_ca)
temp$n <- sapply(pl_ca, function(x,y) nrow(filter(y, Tumor_Sample_Barcode == x)), y = trim_ca)
temp$Cohort <- rep("TCGA Caucasian", length(pl_ca))
temp$ID <- pl_ca
df <- temp

temp <- data.frame(row.names = pl_as)
temp$n <- sapply(pl_as, function(x,y) nrow(filter(y, Tumor_Sample_Barcode == x)), y = trim_as)
temp$Cohort <- rep("TCGA Asian", length(pl_as))
temp$ID <- pl_as
df <- rbind(df, temp)

temp <- data.frame(row.names = pl_ih)
temp$n <- sapply(pl_ih, function(x,y) nrow(filter(y, Tumor_Sample_Barcode == x)), y = trim_ih)
temp$Cohort <- rep("Inhouse Asian", length(pl_ih))
temp$ID <- pl_ih
df <- rbind(df, temp)

temp <- data.frame(row.names = pl_ic)
temp$n <- sapply(pl_ic, function(x,y) nrow(filter(y, Tumor_Sample_Barcode == x)), y = trim_ic)
temp$Cohort <- rep("ICGC", length(pl_ic))
temp$ID <- pl_ic
df <- rbind(df, temp)

temp <- df %>% filter(Cohort != "TCGA Caucasian")
temp$Cohort <- "Combined Asian"
df <- rbind(df, temp)

temp <- df %>% filter(ID %in% p_com_hpv_neg)
temp$Cohort <- "Combined Asian HPV-negative"
df <- rbind(df, temp)

df$Cohort <- factor(df$Cohort, ordered = TRUE, 
                    levels = c("TCGA Caucasian", "TCGA Caucasian HPV positive", "TCGA Caucasian HPV negative",
                               "ICGC", "ICGC HPV positive", "ICGC HPV negative", "TCGA Asian", "Inhouse Asian", "Combined Asian", "Combined Asian HPV-negative"))

#Plot box plot with dot plot
#ggplot(df, aes(x=Cohort, y=n)) + 
#  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +
#  geom_jitter(shape=16, position=position_jitter(0.2))
#p

#theme_set(theme_bw())

#no color
#no significance
####-----####
tiff(file = "test.tiff", height=5, width=15, units="in", res=100, compression="lzw")

#compare_means(n ~ Cohort,  data = df[df$Cohort %in% c("TCGA Caucasian", "Inhouse Asian", "TCGA Asian", "ICGC"),], ref.group = "TCGA Caucasian",
#              method = "t.test")
#my_comparisons <- list( c("TCGA Caucasian", "Inhouse Asian"), c("TCGA Caucasian", "TCGA Asian"), c("TCGA Caucasian", "ICGC") )

p1 <- ggboxplot(df[df$Cohort %in% c("TCGA Caucasian", "Inhouse Asian", "TCGA Asian", "ICGC", "Combined Asian"),], 
             x="Cohort", y="n", color = "Cohort", add = "jitter", palette = "jco") + 
  rotate_x_text(angle = 45)+
  #geom_boxplot() + 
  #geom_dotplot(binaxis='y', 
  #             stackdir='center', 
  #             stackratio=1.5, 
  #             binwidth = 15,
  #             dotsize = 1.3) +
  stat_compare_means(method = "t.test", ref.group = "TCGA Caucasian", label = "p.signif") +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(2000, 3000, 4000)) +
  theme_classic() +
  theme(#axis.text.x = element_text(angle=65, size = 10),
        plot.title = element_text(size = 13, face = "bold", hjust=0.5, vjust = 0),
        legend.position = "none") + 
  labs(title="All tumors", 
       y="Somatic mutation count per tumor",
       x = "")
p2 <- ggplot(df[df$Cohort %in% c("TCGA Caucasian HPV negative", "Inhouse Asian", "TCGA Asian", "ICGC HPV negative", "Combined Asian HPV-negative"),], aes(x=Cohort, y=n)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               stackratio=1.5, 
               binwidth = 15,
               dotsize = 1.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=65, vjust=0.5, size = 10),
        plot.title = element_text(size = 13, face = "bold", hjust=0.5, vjust = 0),) + 
  labs(title="HPV-negative tumors", 
       y="Somatic mutation count per tumor",
       x = "")

p3 <- ggplot(df[df$Cohort %in% c("TCGA Caucasian HPV positive",  "ICGC HPV positive"),], aes(x=Cohort, y=n)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               stackratio=1.5, 
               binwidth = 15,
               dotsize = 1.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=65, vjust=0.5, size = 10),
        plot.title = element_text(size = 13, face = "bold", hjust=0.5, vjust = 0),) + 
  labs(title="HPV-positive tumors", 
       y="Somatic mutation count per tumor",
       x = "")

figure <- ggarrange(p1, p2, p3, align = "h",
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

annotate_figure(figure,
                top = text_grob("Somatic mutation burden", face = "bold", size = 20, vjust = 0.3)
                #bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                #                   hjust = 1, x = 1, face = "italic", size = 10),
                #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                #right = "I'm done, thanks :-)!",
                #fig.lab = "Figure 1", fig.lab.face = "bold"
)

dev.off()
####-----####



#color
#significance
####-----####
tiff(file = "test.tiff", height=5, width=20, units="in", res=100, compression="lzw")

#compare_means(n ~ Cohort,  data = df[df$Cohort %in% c("TCGA Caucasian", "Inhouse Asian", "TCGA Asian", "ICGC"),], ref.group = "TCGA Caucasian",
#              method = "t.test")
#my_comparisons <- list( c("TCGA Caucasian", "Inhouse Asian"), c("TCGA Caucasian", "TCGA Asian"), c("TCGA Caucasian", "ICGC") )

p1 <- ggboxplot(df[df$Cohort %in% c("TCGA Caucasian", "Inhouse Asian", "TCGA Asian", "ICGC","Combined Asian"),], 
                x="Cohort", y="n", color = "Cohort", add = "jitter", palette = "jco") + 
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "t.test", ref.group = "TCGA Caucasian", label = "p.signif") +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(2000, 3000, 4000)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=20, size = 10, vjust = 0.5),
    plot.title = element_text(size = 13, face = "bold", hjust=0.5, vjust = 0),
    legend.position = "none") + 
  labs(title="All tumors", 
       y="Somatic mutation count per tumor",
       x = "")

p2 <- ggboxplot(df[df$Cohort %in% c("TCGA Caucasian HPV negative", "Inhouse Asian", "TCGA Asian", "ICGC HPV negative", "Combined Asian HPV-negative"),], 
                x="Cohort", y="n", color = "Cohort", add = "jitter", palette = "jco") + 
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "t.test", ref.group = "TCGA Caucasian HPV negative", label = "p.signif") +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(2000, 3000, 4000)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=20, size = 10, vjust = 0.5),
    plot.title = element_text(size = 13, face = "bold", hjust=0.5, vjust = 0),
    legend.position = "none") + 
  labs(title="HPV-negative tumors", 
       y="Somatic mutation count per tumor",
       x = "")

p3 <- ggboxplot(df[df$Cohort %in% c("TCGA Caucasian HPV positive",  "ICGC HPV positive"),], 
                x="Cohort", y="n", color = "Cohort", add = "jitter", palette = "jco") + 
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "t.test", ref.group = "TCGA Caucasian HPV positive", label = "p.signif") +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(2000, 3000, 4000)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, size = 10, vjust = 0.5),
    plot.title = element_text(size = 13, face = "bold", hjust=0.5, vjust = 0),
    legend.position = "none") + 
  labs(title="HPV-positive tumors", 
       y="Somatic mutation count per tumor",
       x = "")

figure <- ggarrange(p1, p2, p3, align = "h",
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)

annotate_figure(figure,
                top = text_grob("Somatic mutation burden", face = "bold", size = 20, vjust = 0.3)
                #bottom = text_grob("Data source: \n mtcars data set", color = "blue",
                #                   hjust = 1, x = 1, face = "italic", size = 10),
                #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                #right = "I'm done, thanks :-)!",
                #fig.lab = "Figure 1", fig.lab.face = "bold"
)

dev.off()
####-----####
