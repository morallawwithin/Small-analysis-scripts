library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggbeeswarm)
library(wesanderson)
library(rstatix)
col_wes<-wes_palette("Zissou1",6, type = "continuous")
df<-data.frame("V0.5"=c(wt[,4],V381Y[,4],kcna1[,4],kcna1_kcna2[,4]),
               "k"=c(wt[,5],V381Y[,5],kcna1[,5],kcna1_kcna2[,5]),
               "g"=c(wt[,3],V381Y[,3],kcna1[,3],kcna1_kcna2[,3]),
               "c"=c(wt[,6],V381Y[,6],kcna1[,6],kcna1_kcna2[,6]),
               "group"=factor(c(rep("KCNA2",length(wt[,1])),
                                rep("KCNA2 V381Y",length(V381Y[,1])),
                                rep("KCNA1",length(kcna1[,1])),
                                rep("KCNA1+KCNA2",length(kcna1_kcna2[,1]))
                                )))
for (i in 1:4){
  df[,i]<-as.numeric(df[,i])
}
df$group<-factor(df$group,c("KCNA2","KCNA2 V381Y","KCNA1","KCNA1+KCNA2"))

ggplot(data=df,aes(x=group, y=V0.5, group=group, fill=group,shape = group))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 0.25,  col="black", size=1) +
  geom_beeswarm(cex=5,size=4)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'errorbar', size=1) +
  scale_colour_manual(values = c("#EC7A05","#F11B00","#3A9AB2","#BDC881")) +
  scale_fill_manual(values = c("#EC7A05","#F11B00","#3A9AB2","#BDC881")) +
  scale_shape_manual (values =c(21,22,23,24))+
  theme_prism(base_size = 14)
  

t.test(data=df,V0.5~group)
pairwise_t_test(V0.5 ~ group, p.adjust.method = "bonferroni", data=df)
summary<-df %>%
  group_by( group) %>% 
  summarise(meanV = mean(V0.5),
            sdV = sd(V0.5),
            nV = n(),
            meank = mean(k),
            sdk = sd(k),
            meang = mean(g),
            sdg = sd(g),
            meanc = mean(c),
            sdc = sd(c)
            ) %>%
  mutate(SEMV = sdV/sqrt(nV),
         SEMk = sdk/sqrt(nV),
         SEMg = sdg/sqrt(nV),
         SEMc = sdc/sqrt(nV))

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/activation_2025.png", width = 5, height = 4)


df_act<-rbind(wt_act,V381Y_act,kcna1_act,kcna1_kcna2_act)
df_act$group<-factor(c(rep("KCNA2 wt",length(wt_act$tail_norm)),
                       rep("KCNA2 V381Y",length(V381Y_act$tail_norm)),
                       rep("KCNA1",length(kcna1_act$tail_norm)),
                       rep("KCNA1+KCNA2",length(kcna1_kcna2_act$tail_norm))))
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
  scale_colour_manual(values = c( "#EC7A05","#F11B00","#3A9AB2","#BDC881")) +
  scale_fill_manual(values = c("#EC7A05","#F11B00","#3A9AB2","#BDC881")) +
  scale_shape_manual (values =c(21,22,23,24))+
  xlim(c(-70,60))+
  theme_prism(base_size = 12)+
  xlab("membrane potential [mV]") + ylab("norm. conductance")
ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/activation_curve_all_tails.svg", width = 5, height = 3)
