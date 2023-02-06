library(tidyverse)
library(multimode)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/Caspr1-Kcna2-staining")
df <-
  list.files(path = "D:/Peter/Analysis/KCNA2/P405L_Mice/Caspr1-Kcna2-staining", pattern = "*.csv") %>% 
  map_df(~read_csv(.))
df$position<-rep(1:351,length(df$Line_Intensity)/351)/35.1
df$Line_Intensity<-df$Line_Intensity/max(df$Line_Intensity)
df$staining<-as.factor(df$staining)
ggplot(df_k, aes(position,Line_Intensity, colour=staining))+
  geom_point()+
  geom_smooth(method=lm)
df_k<-filter(df, staining=="anti-kv1.2")
df_c<-filter(df, staining=="anti-caspr")
kcna2<-sapply(unique(df$number), function(x){
df_knr<-filter(df_k,number==x)
df_k0nr<-filter(df_k,Line_Intensity>median(df_knr$Line_Intensity))
locmodes(df_k0nr$position,mod0=2,display=TRUE)})
names(kcna2[,1:length(unique(df$number))])<-unique(df$number)

caspr<-sapply(unique(df$number), function(x){
df_cnr<-filter(df_c,number==x)
df_c0nr<-filter(df_c,Line_Intensity>median(df_cnr$Line_Intensity))
locmodes(df_c0nr$position,mod0=2,display=TRUE)})