library(tidyverse)
library(readxl)
library(ggprism)
library(ggbeeswarm)
library(lsmeans)
library(lme4)
library(ggalluvial)
mainDir<-"D:/Peter/Analysis/KCNA2/P405L_Mice/staining/Golgi" #ordner mit files
subDir<- "0827S1_Golgi" #name des files

dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
data <- read_excel(paste0(mainDir,"/",subDir,".xlsx"))
colnames(data)<-data[7,]
data<-data[-(1:7),]
data<-data[!is.na(data$Mouse),]
data$genotype<-data$Mouse
data$genotype[data$genotype %in% c("0 1","0 4","1 3","0 6","1 5","0 2")]<-"Kcna2+/+"
data$genotype[!(data$genotype=="Kcna2+/+")]<-"Kcna2+/P405L"
data$mouse_neuron<-paste0(data$Mouse,"_",data$`dendr. segment`)

data$Type<-as.factor(data$Type)
data$Length_segment<-as.numeric(data$Length_segment)
data$LWR<-as.numeric(data$LWR)

data_summary<-  data%>%
  group_by(mouse_neuron,genotype,Type,Mouse) %>%
  summarise(
    sd.LWR = sd(LWR, na.rm = TRUE),
    mn.LWR = mean(LWR),
    n_type=n(),
    type_dens=n()/Length_segment,
    .groups="keep"
  )
data_summary<-data_summary[!duplicated(data_summary),]

data_summary$n_type[data_summary$Type=="branched"]<-data_summary$n_type[data_summary$Type=="branched"]/2

data_summary$Type<-factor(data_summary$Type,levels = c("filopodia","long thin","thin","stubby","mushroom","branched"))

P1<-ggplot(data_summary,aes(Type,mn.LWR,color=genotype))+
  geom_boxplot(aes(fill=genotype))+
  geom_point(aes(group=genotype),position=position_dodge(width=0.75))+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  ylab("LWR")+xlab("spine types")+
  theme_prism(base_size = 14)
ggsave(P1,width = 10, height=4,
       file="LWR.png")
P1
lwr<-lmer(mn.LWR~Type*genotype+(1|Mouse),data=data_summary)
summary(lwr)
lsmeans(lwr, pairwise ~ genotype | Type, adjust = "tukey")


P2<-ggplot(data_summary,aes(Type,type_dens,color=genotype))+
  geom_boxplot(aes(fill=genotype))+
  geom_point(aes(group=genotype),position=position_dodge(width=0.75))+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  ylab("spine density [1/µm]")+xlab("spine types")+
  theme_prism(base_size = 14)
P2
type_dens<-lmer(type_dens~Type*genotype+(1|Mouse),data=data_summary)
lsmeans(type_dens, pairwise ~ genotype | Type, adjust = "tukey")
ggsave(P2,width = 10, height=4,
       file="spine_density.png")

data_summary2<-  data %>% 
  group_by(genotype, Type) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)*100) 

data_summary2$Type<-factor(data_summary2$Type,levels = c("filopodia","long thin","thin","stubby","mushroom","branched"))  

P3<-ggplot(data_summary2, aes(genotype,perc, fill = Type)) +
  geom_flow(aes(alluvium = Type), alpha= .2, color = "white",
            curve_type = "linear", 
            width = .5) +
  geom_bar(stat="identity", width = 0.5, color = "white") +
  scale_fill_brewer(palette = "RdBu")+
  ylab("percentage")+
  theme_prism(base_size = 14)+
  theme(legend.position = "none")
ggsave(P3,width = 5, height=4,
       file="spine_composition.svg")  
P3  


data_summary<-  data_summary %>% 
  group_by(mouse_neuron) %>% 
  mutate(perc = n_type/sum(n_type)*100) 


P4<-ggplot(data_summary,aes(Type,perc,color=genotype))+
  geom_boxplot(aes(fill=genotype))+
  geom_point(aes(group=genotype),position=position_dodge(width=0.75))+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  ylab("percentage of spines per neuron [%]")+xlab("spine types")+
  theme_prism(base_size = 14)
ggsave(P4,width = 10, height=4,
       file="spine_percentage.png") 
P4
type_prop<-lmer(perc~Type*genotype+(1|Mouse),data=data_summary)
summary(type_prop)
lsmeans(type_prop, pairwise ~ genotype | Type, adjust = "tukey")
