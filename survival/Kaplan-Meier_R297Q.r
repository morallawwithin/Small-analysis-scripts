library(survival)
library(survminer)
library(dplyr)
library(ggsurvey)
library(svglite)
library(readxl)
library(tidyverse)
library(ggprism)

setwd("G:/Arbeit/Survival")
data<-read_xlsx("R297Q_20241120.xlsx")
data$Genotype<-as.factor(data$Genotype)
data<-data[data$Genotype %in% c("C57Bl/6","SWISS","SWISS + 4 AP"),]
#data<-data[data$Sex=="m",]
surv_object <- Surv(time = data$Lifespan, event = data$Endpoint)
fit1<-survfit(surv_object ~ Genotype, data = data)
survp<-
  ggsurvplot(fit1, data = data, color="strata",
             censor=T,censor.shape=124,censor.size=2,size = 1,
             legend=c(0.8,0.8),legend.title= element_blank(),
             legend.labs=c("C57Bl/6","SWISS","SWISS + 4 AP"),
             xlim = c(0, 150),ylim=c(0,1),xlab="Time (d)",
             ggtheme = theme_prism(base_size = 10))
survp$plot<-survp$plot +   
  scale_x_continuous(expand = c(0, 0, .05, 0)) +
  scale_y_continuous(expand = c(0, 0, .05, 0))+
  theme(legend.spacing.y = unit(-0.2, 'cm'))  +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text.x = element_text(face="plain"))+
  theme(axis.text.y = element_text(face="plain"))+
  scale_color_manual(values = c("lightgrey","red","darkred"))
survp$plot

ggsave("survplot_R297Q_20241120.png", survp$plot, height = 7.5, width = 7.5, units="cm")


survdiff(surv_object ~ Genotype, data = data)
