library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(RColorBrewer)
library(lsmeans)
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

cell01<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_7",
                               pattern="_0(09[1-9]|10[0-1])"),
                           sep=""),
                    condition ), ncol = 2)
cell02<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(4[5-9]|5[0-5])"),
                           sep=""),
                     condition ), ncol = 2)
cell03<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(5[8-9]|6[0-8])"),
                           sep=""),
                     condition ), ncol = 2)
cell04<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(7[0-9]|80)"),
                           sep=""),
                     condition ), ncol = 2)
cell05<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(8[2-9]|9[0-2])"),
                           sep=""),
                     condition ), ncol = 2)
cell06<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12",
                                      pattern="_00(0[1-9]|1[0-1])"),
                           sep=""),
                     condition ), ncol = 2)
cell07<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12",
                                      pattern="_00(3[6-9]|4[0-6])"),
                           sep=""),
                     condition ), ncol = 2)

cellname<-c(paste("cell0",c(1:7),sep=""))
cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07)
cell_values<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(cell_values)<-c("cell","condition","amplitude","v1/2","tau")
cell_values_act<-data.frame(matrix(ncol = 6, nrow = 0))
colnames(cell_values_act) <- c("cell", "condition", "voltage", "cond. norm.", "tail. curr.", "tau")

for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(rep(cellname[i],nrow(curr_cell)),
                       curr_cell[,2],
                       rep(0, nrow(curr_cell)),
                       rep(0, nrow(curr_cell)),
                       rep(0, nrow(curr_cell)))
  
  for (s in 1:nrow(curr_cell)){
    
    sweepnr<-16
    data<-readABF(curr_cell[s,1])
    cond<-rep(0, sweepnr)
    tail_curr<-rep(0, sweepnr)
    curr<-rep(0, sweepnr)
    volt<-10*(-8:7)
    tau <-rep(0, sweepnr)
    for (ii in 1:sweepnr){
      sweep.data<-as.data.frame(data,sweep=ii)
      sweep.max<-sweep.data[c(500:5000),c(1,3)]
      curr[ii]<-max(sweep.max)
      cond[ii]<-max(sweep.max/(volt[i]+95))
      sweep.tail<-sweep.data[c(20000:21000),c(1,3)]
      tail_curr[ii]<-min(sweep.tail)
      # --- Tau of Inactivation ---
      # Fit decay after peak
      peak_idx <- which(sweep.data[,3]==max(sweep.max))
      decay_data <- sweep.data[peak_idx:20000,c(1,3)]
      #decay_data[,2]<-decay_data[,2]-min(decay_data[,2])
      colnames(decay_data)<-c("time","current")
      try({
        fit <- nlsLM(current ~ A * exp(time / tau) + C,
                     start = list(A = decay_data$current[1], tau = -0.5, C = min(decay_data$current)),
                     control = nls.lm.control(maxiter = 500), data=decay_data)
        tau[ii] <- coef(fit)["tau"]
        decay_data$fit<-predict(fit)
        if(ii==12){print(
        ggplot(decay_data, aes(x = time)) +
          geom_line(aes(y = current), color = "blue", size = 1, alpha = 0.6) +
          geom_line(aes(y = fit), color = "red", size = 1) +
          theme_minimal()
        )}
      }, silent = TRUE)
      
      
    }
    cond_norm<-cond/max(cond)
    tail_norm<-tail_curr/min(tail_curr)
    cell_values_act <- rbind(cell_values_act,
                             cbind(rep(cellname[i], sweepnr),
                                   curr_cell[s,2],
                                   volt,
                                   cond_norm,
                                   tail_norm,
                                   tau))
    #cell_values[s,4]<-coef(model)[2]
    cell_values_i[s,3]<-curr[12] #40 mV
    cell_values_i[s,5]<-tau[12] #40 mV
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
cell_values<-cell_values[!cell_values$Condition=="base",]
colnames(cell_values)<-c("cell"  ,      "ident"   ,     "raw.Amp"   ,     "norm.Amp"   ,"Tau",     "Condition", "Time")
cell_values$norm.Amp<-as.numeric(cell_values$norm.Amp)
cell_values$Time<-as.numeric(cell_values$Time)
cell_values$Tau<-as.numeric(cell_values$Tau)
cell_values$Tau[cell_values$Tau==0]<-NA
#write_xlsx(cell_values,"D:/Peter/Analysis/KCNA2/BMK86-P1/V381Y.xlsx")


ggplot(data=cell_values,aes(x=Time, y=norm.Amp, group=Condition, fill=Condition,shape = Condition))+
  coord_cartesian(clip = 'off',ylim=c(0,1.5), xlim = c(0,12))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 0.25,  col="black") +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'path',  size=1,  aes(col=Condition)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'point',  size=4, aes(col=Condition)) +
  scale_colour_manual(values = c("#3c5396", "#b5595f")) +
  scale_fill_manual(values = c("#bac5e3", "#e6a2a4")) +
  scale_shape_manual (values =c(21,22))+
  ylab(expression('norm. K'[V]*'1 current'))+
  theme_prism(base_size = 12)

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2+KCNA1.svg", width = 3.5, height = 2)

cell_values$Time<-as.character(cell_values$Time)
cell_values$Time<-factor(cell_values$Time, levels = c("0","3","6","9","12"))
kcna12_tox<-lm(norm.Amp~Condition*Time, data = cell_values)
anova(kcna12_tox)
lsmeans(kcna12_tox, pairwise ~ Condition | Time, adjust = "tukey")

summary <- cell_values %>%
  group_by( Condition, Time) %>% 
  summarise(meanAmp = mean(norm.Amp),
            sdAmp = sd(norm.Amp),
            nAMp = n()) %>%
  mutate(SEMAmp = sdAmp/sqrt(nAMp))
#########################################
#Inactivation
cell_values$Tau[21]<-NA

ggplot(data=cell_values,aes(x=Time, y=abs(Tau), group=Condition, fill=Condition,shape = Condition))+
  coord_cartesian(clip = 'off',ylim=c(0,4), xlim = c(0,12))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 0.25,  col="black") +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'path',  size=1,  aes(col=Condition)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'point',  size=4, aes(col=Condition)) +
  scale_colour_manual(values = c("#3c5396", "#b5595f")) +
  scale_fill_manual(values = c("#bac5e3", "#e6a2a4")) +
  scale_shape_manual (values =c(21,22))+
  ylab(expression('tau '[inactivation]*'[s]'))+
  theme_prism(base_size = 12)

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2+KCNA1_tau.svg", width = 3.5, height = 2)
kcna12_tau<-lm(Tau~Condition*Time, data = cell_values)
anova(kcna12_tau)
lsmeans(kcna12_tau, pairwise ~ Condition | Time, adjust = "tukey")

summary <- na.omit(cell_values) %>%
  group_by( Condition, Time) %>% 
  summarise(meanTau = mean(Tau),
            sdTau = sd(Tau),
            nTau = n()) %>%
  mutate(SEMTau = sdTau/sqrt(nTau))

#summary(model)
#sweep$PredictionsNLS <- predict(model)
#ggplot(data=sweep, aes(x=voltage,y=conductance))+
#geom_point()+
#geom_line(aes(x=voltage,y=PredictionsNLS))+
#theme_classic()
##################
#Plot examples
###################

data1<-readABF(cell02[1,1])
data1<-as.data.frame(data1,sweep=12)
data2<-readABF(cell02[6,1])
data2<-as.data.frame(data2,sweep=12)
data3<-readABF(cell02[11,1])
data3<-as.data.frame(data3,sweep=12)
ggplot(data1,aes(`Time [s]`,`IN 0C [pA]`))+
  geom_line(color="#252525",size=2)+
  geom_line(data=data2,color="#3c5396",size=2)+
  geom_line(data=data3,color="#b5595f",size=2)+
  ylim(c(-3500,6000))+
  theme_prism(base_size = 14)+
  theme(axis.title.x=element_blank(),        axis.text.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),        axis.text.y=element_blank(),
        legend.position = "none")
ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA12_example.svg", width = 2, height = 2)
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


ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA1+KCNA2_Toxin_act.png", width = 8, height = 6)
