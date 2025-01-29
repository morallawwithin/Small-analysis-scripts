library(survival)
library(survminer)
library(survcomp)
library(dplyr)
library(ggsurvey)
library(svglite)
library(readxl)
library(tidyverse)
library(ggprism)
#setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/Phenomaster")
setwd("G:/Arbeit/Survival")
data<-read_xlsx("survival_P405L_20230612.xlsx")
data<-data[data$Genotype %in% c("C57Bl/6","F1","F2","F3"),]
data$Genotype[data$Genotype %in% c("F1","F2","F3")]<-"Swiss"
data$Genotype<-as.factor(data$Genotype)
surv_object <- Surv(time = data$Lifespan, event = data$Endpoint)
fit1<-survfit(surv_object ~ Genotype, data = data)
survp<-
  ggsurvplot(fit1, data = data, color="strata",
             censor=T,censor.shape=124,censor.size=2,size = 1,
             legend=c(0.8,0.8),legend.title= element_blank(),
             legend.labs=c("C57Bl/6","Swiss"),
             xlim = c(0, 500),ylim=c(0,1),xlab="Time (d)",
             ggtheme = theme_prism(base_size = 10))
survp$plot<-survp$plot +   
  scale_x_continuous(expand = c(0, 0, .05, 0)) +
  scale_y_continuous(expand = c(0, 0, .05, 0))+
  theme(legend.spacing.y = unit(-0.2, 'cm'))  +
  guides(color = guie_legend(byrow = TRUE))+
  theme(axis.text.x = element_text(face="plain"))+
  theme(axis.text.y = element_text(face="plain"))
survp$plot

no.at.risk(surv_object ~ Genotype, data=data, sub.s = "all", t.step=10, t.end=100)

ggsave("survplot_20240416.svg", survp$plot, height = 7.5, width = 7.5, units="cm")


survdiff(surv_object ~ Genotype, data = data)
