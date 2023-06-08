library(readxl)
library(tidyverse)
data<-read_xlsx("G:/Arbeit/EEG/Auswertung_20230607.xlsx")
data$animal<-as.factor(data$animal)
data$Area<-as.factor(data$Area)
data$Genotype<-factor(data$Genotype,levels = c("wt","P405L"))
data<-data[!is.na(data$`%`),]
data<-data[!(data$Area %in% c("PPAA","TH")),]
areas<-aov(data$`%`~Area*Genotype+Error(animal),data=data)
summary(areas)

ggplot(data, aes(Area,data$`%`,color=Genotype))+
  geom_boxplot(aes(fill=Genotype))+
  geom_point(aes(group=Genotype),position=position_dodge(width=0.75))+
  theme_classic()+
  scale_color_manual(values = c("black","blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  ylab("c-fos positive neurons [%]")
