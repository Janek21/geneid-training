library(ggplot2)

cmp_data<-read.csv("../results/gffcmp/specie_stats.csv", sep="\t")
colnames(cmp_data)
rownames(cmp_data)<-cmp_data$Species

cmp_data<-cmp_data[cmp_data$Species %in% c("Drosophila_melanogaster", "Drosophila_willistoni", "Plasmodium_falciparum", "Plasmodium_vivax"),]
#total precision/sensitivity by species
target_cols<-c("Base.level.Sensitivity", "Base.level.Precision", "Exon.level.Sensitivity", "Exon.level.Precision", "Intron.level.Sensitivity", "Intron.level.Precision", "Transcript.level.Sensitivity", "Transcript.level.Precision")
precSens_long<-reshape(cmp_data, 
                       varying=target_cols,
                       v.names="Value",
                       timevar="Metric",
                       times=target_cols,
                       direction="long")

#separate metric(sensitivity, precision) from level(exon, intron, base)
precSens_long$Level <- gsub("\\.level\\.(Sensitivity|Precision)", "", precSens_long$Metric)
precSens_long$Metric <- gsub(".*\\.level\\.", "", precSens_long$Metric)

#plot
ggplot(precSens_long, aes(x=Species, y=Value, fill=Metric))+geom_col(position="dodge")+
  labs(title="Sensitivity/Precision scoring", y="Percentage score")+
  scale_fill_brewer(palette="Set1")+
  facet_wrap(~Level)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))






#novel/missed
target_cols<-c("Novel.exons....", "Missed.loci....", "Missed.introns....", "Novel.loci....", "Missed.exons....", "Novel.introns....")
mn_long<-reshape(cmp_data, 
                       varying=target_cols,
                       v.names="Value",
                       timevar="Metric",
                       times=target_cols,
                       direction="long")

#separate metric(novel, missing) from level(exon, intron, base)
mn_long$Level <- gsub("Missed\\.|Novel\\.|\\.\\.\\.\\.", "", mn_long$Metric)
mn_long$Metric <- ifelse(grepl("Missed", mn_long$Metric), "Missed", "Novel")

#plot
ggplot(mn_long, aes(x=Species, y=Value, fill=Metric))+geom_col(position="dodge")+
  labs(title="Novel/Missed scoring", y="Percentage score")+
  scale_fill_manual(values = c("Missed" = "#7f8c8d", "Novel" = "#f39c12"))+
  facet_wrap(~Level)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))

