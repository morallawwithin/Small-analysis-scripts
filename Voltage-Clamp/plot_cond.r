library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggbeeswarm)

df<-data.frame("V0.5"=c(wt,V381Y),
               "group"=factor(c(rep("wt",15),
                                rep("V381Y",8))))
df$group<-factor(df$group,c("wt","V381Y"))

ggplot(data=df,aes(x=group, y=V0.5, group=group, fill=group,shape = group))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 0.25,  col="black", size=1) +
  geom_beeswarm(size=4)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'errorbar', size=1) +
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("black", "blue")) +
  scale_shape_manual (values =c(21,22))+
  theme_prism(base_size = 14)
  

t.test(data=df,V0.5~group)

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/activation.png", width = 5, height = 4)
