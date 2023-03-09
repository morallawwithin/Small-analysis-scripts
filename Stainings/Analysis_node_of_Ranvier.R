library(tidyverse)
library(multimode)
library(lme4)
library(lmerTest)
library(cowplot)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/Caspr1-Kcna2-staining")
df <-
  list.files(path = "D:/Peter/Analysis/KCNA2/P405L_Mice/Caspr1-Kcna2-staining", pattern = "*.csv") %>% 
  map_df(~read_csv(.))
df$position<-rep(1:351,length(df$Line_Intensity)/351)/35.1
df$Genotype<-as.factor(df$Genotype)
df$staining<-as.factor(df$staining)
df<-filter(df,Line_Intensity>0)

df_k<-filter(df, staining=="anti-kv1.2")

o<-stack(sapply(unique(df_k$imagename),function(x){
  (df_k$Line_Intensity[df_k$imagename==x]-mean(df_k$Line_Intensity[df_k$imagename==x]))/sd(df_k$Line_Intensity[df_k$imagename==x])
}))
df_k$Line_Intensity<-o$values



####Playing around, don't use
df_kz<-filter(df_k,Line_Intensity>0)
df_kz$position<-as.factor(df_kz$position)
ggplot(filter(df_kz,df_kz$Genotype=="Kcna2 +/+"), aes(position,Line_Intensity))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = quantile(df_kz$Line_Intensity, c(0.1, 0.9)))

  geom_boxplot(data = filter(df_kz,df_kz$Genotype=="Kcna2 +/+"))+
  geom_boxplot(data = filter(df_kz,df_kz$Genotype=="Kcna2 +/P405L"),colour=blue)
k<-lmer(Line_Intensity~position+(1|number)+(1|imagename), data=df_k)
summary(k)
#####

kcna2<-lapply(unique(df$number), function(x){
df_knr<-filter(df_k,number==x)
df_k0nr<-filter(df_knr,Line_Intensity>0)
locmodes(df_k0nr$position,mod0=2,display=TRUE)})
kcna2_v<-matrix(, nrow = length(unique(df$number)), ncol = 2)
for (y in c(1:length(unique(df$number)))) {
  kcna2_v[y,1] <- 5-kcna2[[y]][["locations"]][1]
  kcna2_v[y,2] <- kcna2[[y]][["locations"]][3]-5 
}
kcna2_v

df_k$nodeid<-paste(df_k$imagename,df_k$pointsno, sep = "_")

nod<-sapply(unique(df_k$nodeid), function(x){
  df_kn<-filter(df_k,nodeid==x)
  sum(as.numeric(any(df_kn$Line_Intensity>1.5&df_kn$position<5)),as.numeric(any(df_kn$Line_Intensity>1.5&df_kn$position>5)))
})
nodes<-data.frame("animal"=as.factor(df_k$number[match(unique(df_k$nodeid),df_k$nodeid)]), 
                  "Genotype"=as.factor(df_k$Genotype[match(unique(df_k$nodeid),df_k$nodeid)]),
                  "nodes_nr"=nod)
a<-ggplot(nodes,aes(Genotype,nodes_nr,color=Genotype,fill=Genotype))+
  geom_violin()+
  geom_jitter(height=0)+
  scale_color_manual(values = c("black","blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  ylab("Kv1.2 patches per node")+
  theme_half_open()+
  theme(legend.position = "none")

ggsave(a, filename = "kcna2_patches_per_node.svg", height = 4, width=2)

wilcox.test(nodes_nr~Genotype, data = nodes, correct = FALSE)
node_model<-lmer(nodes_nr~Genotype+(1|animal), data=nodes)
summary(node_model)
#plot
df_kn<-filter(df_k,nodeid=="P405L_366_3.lif - 63x3_1")
ggplot(df_kn, aes(position,Line_Intensity))+
  geom_point()+
  ylim(-2, 6)+
  xlab("Position [µm]")+
  ylab("Intensity [z-score]")
#figure examples
#2 nodes
df_kn<-filter(df_k,nodeid=="P405L_366_3.lif - 63x2_7")
ggsave(filename ="P405L_366_3.lif - 63x2_7.svg", height = 1.5, width=4)
#1 node
df_kn<-filter(df_k,nodeid=="P405L_366_3.lif - 63x2_2")
ggsave(filename ="P405L_366_3.lif - 63x2_2.svg", height = 1.5, width=4)
#0 node
df_kn<-filter(df_k,nodeid=="P405L_366_3.lif - 63x2_6")
ggsave(filename ="P405L_366_3.lif - 63x2_5.svg", height = 1.5, width=4)
#length
len<-t(sapply(unique(df_k$nodeid), function(x){
  df_kn<-filter(df_k,nodeid==x)
  c(sum(df_kn$Line_Intensity>1&df_kn$position<5,na.rm=TRUE),sum(df_kn$Line_Intensity>1&df_kn$position>5,na.rm=TRUE))

}))
len<-c(len[,1],len[,2])/35
nodelen<-data.frame("animal"=rep(as.factor(df_k$number[match(unique(df_k$nodeid),df_k$nodeid)]),2), 
                  "Genotype"=rep(as.factor(df_k$Genotype[match(unique(df_k$nodeid),df_k$nodeid)]),2),
                  "nodes_len"=len)
ggplot(nodelen,aes(Genotype,nodes_len,color=Genotype))+
  geom_violin()+
  geom_jitter(height=0)+
  scale_color_manual(values = c("black","blue"))+
  ylab("juxtaparanode length [µm]")
wilcox.test(nodes_len~Genotype, data = nodelen)

##################################################################################################################
df_c<-filter(df, staining=="anti-caspr")
o<-stack(sapply(unique(df_c$imagename),function(x){
  (df_c$Line_Intensity[df_c$imagename==x]-mean(df_c$Line_Intensity[df_c$imagename==x]))/sd(df_c$Line_Intensity[df_c$imagename==x])
}))
df_c$Line_Intensity<-o$values
caspr<-lapply(unique(df$number), function(x){
df_cnr<-filter(df_c,number==x)
df_c0nr<-filter(df_cnr,Line_Intensity>0)
locmodes(df_c0nr$position,mod0=2,display=TRUE)})
caspr_v<-matrix(, nrow = length(unique(df$number)), ncol = 2)
for (y in c(1:length(unique(df$number)))) {
  caspr_v[y,1] <- 5-caspr[[y]][["locations"]][1]
  caspr_v[y,2] <- caspr[[y]][["locations"]][3]-5 
}
caspr_v
