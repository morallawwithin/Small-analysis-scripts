library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(RColorBrewer)
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

cell01<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3",
                               pattern="_00(1[0-9]|20)"),
                           sep=""),
                    condition ), ncol = 2)

cell02<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3",
                               pattern="_00(2[6-9]|3[0-6])"),
                           sep=""),
                    condition ), ncol = 2)

cell03<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3",
                               pattern="_00(4[2-9]|5[0-2])"),
                           sep=""),
                    condition ), ncol = 2)

cell04<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_3",
                               pattern="_00(5[4-9]|6[0-4])"),
                           sep=""),
                    condition ), ncol = 2)
cell05<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4",
                               pattern="_00(09|1[0-9])"),
                     sep=""),
                    condition ), ncol = 2)
cell06<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4",
                               pattern="_00(3[2-9]|4[0-2])"),
                           sep=""),
                    condition ), ncol = 2)
cell07<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4",
                               pattern="_00(49|5[0-9])"),
                           sep=""),
                    condition ), ncol = 2)
cell08<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_4",
                               pattern="_00(6[3-9]|7[0-3])"),
                           sep=""),
                    condition ), ncol = 2)
cell09<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_5/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_5",
                               pattern="_00(0[4-9]|1[0-4])"),
                           sep=""),
                    condition ), ncol = 2)
cell10<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_8/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_8",
                               pattern="_00(4[7-9]|5[0-7])"),
                     sep=""),
                    condition ), ncol = 2)
cell11<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_8/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_8",
                               pattern="_00(6[4-9]|7[0-4])"),
                           sep=""),
                    condition ), ncol = 2)

cell12<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_8/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_8",
                               pattern="_00(8[6-9]|9[0-6])"),
                           sep=""),
                    condition ), ncol = 2)
cell13<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11",
                               pattern="_00(0[5-9]|1[0-5])"),
                           sep=""),
                    condition ), ncol = 2)
cell14<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11",
                               pattern="_00(2[0-9]|30)"),
                           sep=""),
                    condition ), ncol = 2)
cell15<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11",
                               pattern="_01(1[4-9]|2[0-4])"),
                           sep=""),
                    condition ), ncol = 2)

cellname<-c(paste("cell0",c(1:9),sep=""),paste("cell",c(10:15),sep=""))
cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, cell11, cell12, cell13, cell14, cell15)
cell_values<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(cell_values)<-c("cell","condition","amplitude","v1/2")
cell_values_act<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(cell_values_act)<-c("cell","condition","voltage","cond. norm.","tail. curr.")

for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(rep(cellname[i],nrow(curr_cell)),
                       curr_cell[,2],
                       rep(0, nrow(curr_cell)),
                       rep(0, nrow(curr_cell)))
  
  for (s in 1:nrow(curr_cell)){
    
    sweepnr<-16
    data<-readABF(curr_cell[s,1])
    cond<-rep(0, sweepnr)
    tail_curr<-rep(0, sweepnr)
    curr<-rep(0, sweepnr)
    volt<-10*(-8:7)
    
    for (ii in 1:sweepnr){
      sweep.data<-as.data.frame(data,sweep=ii)
      sweep.max<-sweep.data[c(500:5000),c(1,3)]
      curr[ii]<-max(sweep.max)
      cond[ii]<-max(sweep.max/(volt[i]+95))
      sweep.tail<-sweep.data[c(20000:21000),c(1,3)]
      tail_curr[ii]<-min(sweep.tail)
    }
    cond_norm<-cond/max(cond)
    tail_norm<-tail_curr/min(tail_curr)
    cell_values_act<-rbind(cell_values_act,
                           cbind(rep(cellname[i],sweepnr),
                                 curr_cell[s,2],
                                 volt,
                                 cond_norm,
                                 tail_norm)) 
    #cell_values[s,4]<-coef(model)[2]
    cell_values_i[s,3]<-curr[12] #40 mV
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

#write_xlsx(cell_values,"D:/Peter/Analysis/KCNA2/BMK86-P1/V381Y.xlsx")


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

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2.png", width = 5, height = 3)

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
model_extra <- nls(cond_norm ~ activation(myg,myVhalf,myk,myc,volt), data=filter(df_act,CondTime=="base"), start=list(myg=1,myVhalf=-12,myk=13,myc=0),control = nls.control(maxiter = 400))
model_tox<- nls(cond_norm ~ activation(myg,myVhalf,myk,myc,volt), data=filter(df_act,CondTime=="toxin.12"), start=list(myg=1,myVhalf=-12,myk=12,myc=0),control = nls.control(maxiter = 400))

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
                                 start=list(myg=1,myVhalf=-12,myk=13,myc=0)), 
              data = df_act,
              se = FALSE,
  )+
  
  scale_colour_manual(values = c("black", mypal(5),mypal2(5))) +
  scale_fill_manual(values = c("black", mypal(5),mypal2(5))) +
  #scale_shape_manual (values =c(21,22,23,24))+
  #xlim(c(-60,70))+
  theme_prism(base_size = 14)+
  xlab("memb. pot. [mV]") + ylab("norm. cond.")

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2_Toxin_act.png", width = 8, height = 6)
