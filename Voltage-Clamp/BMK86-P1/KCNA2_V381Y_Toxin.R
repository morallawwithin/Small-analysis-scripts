library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(RColorBrewer)
library(lsmeans)
cell1<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0003.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0004.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0005.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0006.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0007.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0008.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0009.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0010.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0011.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0012.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0013.abf","toxin.12"),
                nrow=2)

cell2<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0013a.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0014.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0015.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0016.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0017.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0018.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0019.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0020.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0021.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0022.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0023.abf","toxin.12"),
              nrow=2)

cell3<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0000.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0001.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0002.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0003.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0004.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0005.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0006.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0007.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0008.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0009.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0010.abf","toxin.12"),
              nrow=2)

cell4<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0007.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0008.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0009.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0010.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0011.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0012.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0013.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0014.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0015.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0016.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0017.abf","toxin.12"),
              nrow=2)

cell5<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0007.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0008.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0009.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0010.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0011.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0012.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0013.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0014.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0015.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0016.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0017.abf","toxin.12"),
              nrow=2)

cell6<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0003_#1_Baseline.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0004_#1_WOBMK86-p1.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0005_#1_WOBMK86-p1.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0006_#1_WOBMK86-p1.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0007_#1_WOBMK86-p1.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0008#1_1_BMK86-p1.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0009_#1_2_BMK86-p1.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0010_#1_3_BMK86-p1.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0011_#1_4_BMK86-p1.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0012_#1_5_BMK86-p1.abf","toxin.12"),
              nrow=2)

cell7<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0014_#2_Baseline.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0015_#2_WOBMK86-p1.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0016_#2_WOBMK86-p1.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0017_#2_WOBMK86-p1.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0018_#2_WOBMK86-p1.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0019_#2_BMK86-p1.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0020_#2_BMK86-p1.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0021_#2_BMK86-p1.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0022_#2_BMK86-p1.abf","toxin.9"),
              nrow=2)

cell8<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0012.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0013.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0014.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0015.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0016.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0017.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0018.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0019.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0020.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0021.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0022.abf","toxin.12"),
              nrow=2)
#curr_cell<-cell5
cell_values<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(cell_values)<-c("cell","condition","amplitude","v1/2")
cell_values_act<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(cell_values_act)<-c("cell","condition","voltage","cond. norm.","tail. curr.")
cells<-list(cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8)
cellname<-c("cell1","cell2","cell3","cell4","cell5","cell6","cell7","cell8")
for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(rep(cellname[i],ncol(curr_cell)),
                       curr_cell[2,],
                       rep(0, ncol(curr_cell)),
                       rep(0, ncol(curr_cell)))

  for (s in 1:ncol(curr_cell)){

    data<-readABF(curr_cell[1,s])
    cond<-rep(0, 11)
    tail_curr<-rep(0, 11)
    curr<-rep(0, 11)
    volt<-(-6:4)*10

    for (t in 1:11){
      sweep.data<-as.data.frame(data,sweep=t)
      sweep.max<-sweep.data[c(780:45000),c(1,4)]
      curr[t]<-max(sweep.max)
      cond[t]<-max(sweep.max/(volt[t]+95))
      sweep.tail<-sweep.data[c(45000:46000),c(1,4)]
     tail_curr[t]<-min(sweep.tail)
      }
  cond_norm<-cond/max(cond)
  tail_norm<-tail_curr/min(tail_curr)
  cell_values_act<-rbind(cell_values_act,
                         cbind(rep(cellname[i],11),
                               curr_cell[2,s],
                               volt,
                               cond_norm,
                               tail_norm)) 
  #cell_values[s,4]<-coef(model)[2]
  cell_values_i[s,3]<-curr[10]
  }
  
cell_values<-rbind(cell_values,cell_values_i)
}


cell_values[,3]<-as.numeric(cell_values[,3])

for ( i in cellname){

  cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"ext")),4]<-cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"ext")),3]/cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"bas")),3]
  if (any((cell_values[,1]==i)&(str_detect(cell_values[,2],"extra.12")))){
    cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"tox")),4]<-cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"tox")),3]/cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra.12")),3]
  } else{
    cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"tox")),4]<-cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"tox")),3]/cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra.9")),3]
  }
  }

cell_values [c('Condition', 'Time')]<- str_split_fixed(cell_values$V2,"\\.",2)
cell_values$Time<-factor(cell_values$Time, levels = c("0","3","6","9","12"))
cell_values<-cell_values[!cell_values$Condition=="base",]
colnames(cell_values)<-c("cell"  ,      "ident"   ,     "raw.Amp"   ,     "norm.Amp"   ,     "Condition", "Time")
cell_values$norm.Amp<-as.numeric(cell_values$norm.Amp)

ggplot(data=cell_values,aes(x=Time, y=norm.Amp, group=Condition, fill=Condition,shape = Condition))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 0.25,  position = position_dodge(width = 0.2),col="black") +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'path',  size=1, position = position_dodge(width = 0.2), aes(col=Condition)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'point',  size=4, position = position_dodge(width = 0.2)) +
  scale_colour_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  scale_shape_manual (values =c(21,22))+
  ylim(c(0,1.5))+
  theme_prism(base_size = 14)

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/V381Y.png", width = 5, height = 3)

summary <- cell_values %>%
  group_by( Condition, Time) %>% 
  summarise(meanAmp = mean(norm.Amp),
            sdAmp = sd(norm.Amp),
            nAMp = n()) %>%
  mutate(SEMAmp = sdAmp/sqrt(nAMp))

V381Y_tox<-lm(norm.Amp~Condition*Time, data = cell_values)
anova(V381Y_tox)
lsmeans(V381Y_tox, pairwise ~ Condition | Time, adjust = "tukey")

#summary(model)
#sweep$PredictionsNLS <- predict(model)
#ggplot(data=sweep, aes(x=voltage,y=conductance))+
  #geom_point()+
  #geom_line(aes(x=voltage,y=PredictionsNLS))+
  #theme_classic()

################################################################################################
##I/V curve over time 
################################################################################################
condition<-c(
  "base",
  "extra.0",
  "extra.3",
  "extra.6",
  "extra.9",
  "extra.12",
  "toxin.0",
  "toxin.3",
  "toxin.6",
  "toxin.9",
  "toxin.12")
df_act<-cell_values_act
df_act [c('Condition', 'Time')]<- str_split_fixed(df_act$V2,"\\.",2)
names(df_act)[names(df_act) == 'V2'] <- "CondTime"
#df_act<-df_act[!df_act$Condition=="base",]
df_act$Time<-factor(df_act$Time, levels = c("0","3","6","9","12"))
df_act$CondTime<-factor(df_act$CondTime, levels = condition)
df_act$volt<-as.numeric(df_act$volt)
df_act$cond_norm<-as.numeric(df_act$cond_norm)
df_act$tail_norm<-as.numeric(df_act$tail_norm)


activation <- function(g, Vhalf, k,c,V) (g/(1+exp((V-Vhalf)/k))+c)
model_extra <- nls(cond_norm ~ activation(myg,myVhalf,myk,myc,volt), data=filter(df_act,CondTime=="extra.12"), start=list(myg=1,myVhalf=-10,myk=12,myc=0),control = nls.control(maxiter = 400))
model_tox<- nls(cond_norm ~ activation(myg,myVhalf,myk,myc,volt), data=filter(df_act,CondTime=="toxin.12"), start=list(myg=1,myVhalf=-15,myk=12,myc=0),control = nls.control(maxiter = 400))

mypal <- colorRampPalette(brewer.pal(3, "Blues"),bias = 5)
mypal2 <- colorRampPalette(brewer.pal(3, "YlOrRd"),bias = 5)

ggplot(data=df_act,aes(x=volt, y=cond_norm, group=CondTime, fill=CondTime,shape = Condition,col=CondTime))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 2,  size=1) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'point', size=4) +
  #geom_beeswarm(size=1)+
  geom_smooth(method = "nls", 
              method.args = list(formula = y ~ activation(myg,myVhalf,myk,myc,x),
                                 start=list(myg=1,myVhalf=-12,myk=12,myc=0)), 
              data = df_act,
              se = FALSE,
              )+

  scale_colour_manual(values = c("black", mypal(5),mypal2(5))) +
  scale_fill_manual(values = c("black", mypal(5),mypal2(5))) +
  #scale_shape_manual (values =c(21,22,23,24))+
  #xlim(c(-60,70))+
  theme_prism(base_size = 14)+
  xlab("memb. pot. [mV]") + ylab("norm. cond.")

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2_V381Y_Toxin_act.png", width = 8, height = 6)
