library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggbeeswarm)

df<-data.frame("V0.5"=c(wt,V381Y),
               "group"=factor(c(rep("wt",15),
                                rep("V381Y",13))))
df$group<-factor(df$group,c("wt","V381Y"))

ggplot(data=df,aes(x=group, y=V0.5, group=group, fill=group,shape = group))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 0.25,  col="black", size=1) +
  geom_beeswarm(cex=5,size=4)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'errorbar', size=1) +
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("black", "blue")) +
  scale_shape_manual (values =c(21,22))+
  theme_prism(base_size = 14)
  

t.test(data=df,V0.5~group)
summary<-df %>%
  group_by( group) %>% 
  summarise(meanV = mean(V0.5),
            sdV = sd(V0.5),
            nV = n()) %>%
  mutate(SEMV = sdV/sqrt(nV))

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/activation_2025.png", width = 5, height = 4)


df_act<-rbind(wt_act,V381Y_act,kcna1_act,kcna1_kcna2_act)
df_act$group<-factor(c(rep("KCNA2 wt",15*16),
                       rep("KCNA2 V381Y",(7*16+6*11)),
                       rep("KCNA1",11*16),
                       rep("KCNA1+KCNA2",7*16)))
df_act$group<-factor(df_act$group,c("KCNA2 wt","KCNA2 V381Y","KCNA1","KCNA1+KCNA2"))
df_act$volt<-as.numeric(df_act$volt)
df_act$cond_norm<-as.numeric(df_act$cond_norm)
df_act$tail_norm<-as.numeric(df_act$tail_norm)

activation <- function(g, Vhalf, k,c,V) (g/(1+exp((V-Vhalf)/k))+c)
model_kcna2 <- nls(cond_norm ~ activation(myg,myVhalf,myk,myc,volt), data=filter(df_act, group=="wt"), start=list(myg=1,myVhalf=10,myk=12,myc=0),control = nls.control(maxiter = 400))
model_V381Y<- nls(cond_norm ~ activation(myg,myVhalf,myk,myc,volt), data=filter(df_act, group=="V381Y"), start=list(myg=1,myVhalf=10,myk=12,myc=0),control = nls.control(maxiter = 400))



ggplot(data=df_act,aes(x=volt, y=cond_norm, group=group, fill=group,shape = group))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 4,  col="black", size=1) +
  geom_smooth(method = "nls", 
              method.args = list(formula = y ~ activation(myg,myVhalf,myk,myc,x),
                                 start=list(myg=1,myVhalf=0,myk=9,myc=0)), 
              data = df_act,
              se = FALSE,
              aes(color = group))+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'point', size=4) +
  scale_colour_manual(values = c("black", "blue","yellow","green")) +
  scale_fill_manual(values = c("black", "blue","yellow","green")) +
  scale_shape_manual (values =c(21,22,23,24))+
  xlim(c(-80,80))+
  theme_prism(base_size = 14)+
  xlab("memb. pot. [mV]") + ylab("norm. cond.")
ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/activation_curve_all_tails.png", width = 5, height = 4)
