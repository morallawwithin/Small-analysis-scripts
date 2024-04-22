library(readxl)
library(tidyverse)
library(plyr)
library(lsmeans)
library(lme4)
data<-read_xlsx("D:/Peter/Analysis/KCNA2/P405L_Mice/c-fos/Auswertung_final_20230919_PM.xlsx")#"G:/Arbeit/EEG/Auswertung_final_20230919_PM.xlsx")
colnames(data)[6] ="cfos"
data$animal<-sub("#","",data$animal)
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

data_eeg<-read_xlsx("D:/Peter/Analysis/KCNA2/P405L_Mice/EEG/EEG-Data.xlsx")
data_eeg$`Animal-Nr.`<-as.character(data_eeg$`Animal-Nr.`)
for (i in 1:nrow(data)){
  if(any(data$animal[i]==data_eeg$`Animal-Nr.`) ){
   data$seizure_freq[i]<-data_eeg$`Seizure/day`[data_eeg$`Animal-Nr.`==data$animal[i]]
  }else{
    data$seizure_freq[i]<- NA
  }
}

seizure_cfos<-lm(cfos~seizure_freq,data= filter(data,Genotype=="P405L"))
summary(seizure_cfos)


ggplot(data, aes(seizure_freq,cfos,color=Genotype))+
  geom_point(aes(group=Genotype))+
  stat_smooth(method = "lm", col = rgb(191/255,191/255,1,1),data = filter(data,Genotype=="P405L"))+
  theme_classic()+
  scale_color_manual(values = c("black","blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  ylab("c-fos positive neurons [%]")+
  xlab("seizures frequency [1/d]")


