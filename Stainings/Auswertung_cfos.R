library(readxl)
library(tidyverse)
library(plyr)
library(lsmeans)
library(lme4)
data<-read_xlsx("G:/Arbeit/EEG/Auswertung_final_20230919_PM.xlsx")
colnames(data)[6] ="cfos"
data$animal<-as.factor(data$animal)
data$Area<-as.factor(data$Area)
data$Genotype<-factor(data$Genotype, levels=c("WT","P405L"))
data$Genotype<-revalue(data$Genotype,c("WT"="wt","P405L"="P405L"))
data<-data[!is.na(data$cfos),]
data<-data[!(data$Area %in% c("PPAA")),]
areas<-lmer(cfos~Area*Genotype+(1|animal),data=data)
summary(areas)

lsmeans(areas, pairwise ~ Genotype | Area, adjust = "tukey")


ggplot(data, aes(Area,cfos,color=Genotype))+
  geom_boxplot(aes(fill=Genotype))+
  geom_point(aes(group=Genotype),position=position_dodge(width=0.75))+
  theme_classic()+
  scale_color_manual(values = c("black","blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  ylab("c-fos positive neurons [%]")
