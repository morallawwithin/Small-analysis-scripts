library(tidyverse)
library(data.table)
library(ggprism)
library(ggbeeswarm)
library(nlme)
library(lme4)
library(lmerTest)

setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/Cortex_L2&3_PN/EPSC_Spike")

filelist = list.files(pattern = ".*.txt")

data <- rbindlist(sapply(filelist, fread, simplify = FALSE),
                  use.names = TRUE, idcol = "FileName")

data$Genotype <- sapply(X = strsplit(data$FileName, split = "_"), FUN = "[", 3)
data$age<- sapply(X = strsplit(data$FileName, split = "_"), FUN = "[", 4)
data$animal<- sapply(X = strsplit(data$FileName, split = "_"), FUN = "[", 1)
names(data)[names(data) == "0 Measure 1"] <- 'amplitude'
data$Genotype<-factor(data$Genotype, levels = c("wt","P405L"))
data$IF<-unlist(sapply(unique(data$FileName),FUN = function(x){
  1/diff(c(0,data$Time[data$FileName==x]))
}))
########
#data summary statistics
#######
data_summary<-group_by(data, Genotype, age,FileName) %>% 
  summarise(count = n(),
            time = max(Time),
            mnIF=mean(IF),
            mnAmp=mean(amplitude))%>% 
  mutate(freq = (count/time))
data_summary$Genotype<-factor(data_summary$Genotype,levels = c("wt","P405L"))

###############
#amplitude
##############
p1<-ggplot(data,aes(x=abs(amplitude),color=Genotype))+
  geom_step(stat="ecdf", size=1.5)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,100))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1))+
  ylab("cumulative fraction")+
  xlab("absolute amplitude [pA]")+
  theme_prism(base_size = 14)+
  theme(legend.position = "none")  
p1
ggsave(p1,width = 5, height = 4,
       file="EPSC_amp_cum.png")

ggplot(data,aes(Genotype,amplitude,fill=Genotype, col=Genotype))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.05, outliers = FALSE)+
  #geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  #scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
  ylab("amplitude")+
  theme_prism(base_size = 14)+
  theme(legend.position = "none")

p3<-ggplot(data_summary,aes(Genotype,mnAmp,fill=Genotype, col=Genotype))+
  geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  #scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
  ylab("amplitude")+
  theme_prism(base_size = 14)+
  theme(legend.position = "none")  
ggsave(p3,width = 3, height = 4,
       file="EPSC_amp_boxplot.png")


ampl <-lmer(amplitude~Genotype+(1|animal),data=data)
summary(ampl)

##############
#frequency
#############




wilcox.test(mnAmp~Genotype,data_summary)
wilcox.test(mnIF~Genotype,data_summary)

p2<-ggplot(data,aes(x=IF,color=Genotype))+
  geom_step(stat="ecdf", size=1.5)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,200))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1))+
  ylab("cumulative fraction")+
  xlab("instantaneous frequency")+
  theme_prism(base_size = 14,base_family = "Calibri")+
  theme(legend.position = "none")  
p2
ggsave(p2,width = 5, height = 4,
       file="EPSC_freq_cum.png")
ggsave(p2,width = 4, height = 4,
       file="EPSC_freq_cum.svg")

ggplot(data,aes(FileName,IF,fill=Genotype, col=Genotype))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.05, outliers = FALSE)+
  #geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  #scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
  ylab("amplitude")+
  theme_prism(base_size = 14)+
  theme(legend.position = "none")

p4<-ggplot(data_summary,aes(Genotype,mnIF,fill=Genotype, col=Genotype))+
  geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  #scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
  ylab("instantaneous frequency")+
  theme_prism(base_size = 14)+
  theme(legend.position = "none")  
ggsave(p4,width = 3, height = 4,
       file="EPSC_freq_boxplot.png")

p5<-ggplot(data_summary,aes(Genotype,freq,fill=Genotype, col=Genotype))+
  geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  #scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
  ylab("frequency")+
  theme_prism(base_size = 14)+
  theme(legend.position = "none")  
ggsave(p5,width = 3, height = 4,
       file="EPSC_freq_one_value_boxplot.png")
freq2<-lmer(IF~Genotype+(1|FileName),data=data)
summary(freq2)
freq3<-lm(IF~Genotype,data=data)
summary(freq3)
anova(freq3)
