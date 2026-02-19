library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(RColorBrewer)
library(lsmeans)
condition<-c(
  "base",
  "extra1.0",
  "extra1.3",
  "extra1.6",
  "extra1.9",
  "extra1.12",
  "extra2.0",
  "extra2.3",
  "extra2.6",
  "extra2.9",
  "extra2.12")

cell2<-matrix(c("D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711003.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711004.abf","extra1.0",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711005.abf","extra1.3",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711006.abf","extra1.6",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711007.abf","extra1.9",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711009.abf","extra2.3",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711010.abf","extra2.6",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711011.abf","extra2.9",
                "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250711/Kv1.2_WT/25711012.abf","extra2.12"),
              nrow=2)

cell3<-matrix(rbind(c(
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804000.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804001.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804002.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804003.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804004.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804005.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804006.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804007.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804008.abf"
),c(
  "base",
  "extra1.0",  "extra1.3",  "extra1.6",  "extra1.9",  "extra1.12",
  "extra2.0",  "extra2.3",  "extra2.6")),
nrow=2)

cell4<-matrix(rbind(c(

  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804009.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804010.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804011.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804012.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804013.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804014.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804015.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804016.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804017.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804018.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250804/WT/25804019.abf"),
  c(
  "base",
  "extra1.0",  "extra1.3",  "extra1.6",  "extra1.9",  "extra1.12",
  "extra2.0",  "extra2.3",  "extra2.6","extra2.9","extra2.12")),
  nrow=2)

cell5<-matrix(rbind(c(
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805001.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805002.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805003.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805004.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805005.abf",

  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805007.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805008.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805009.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805010.abf",
  "D:/Peter/Data/KCNA2/BMK86/Remy_BMK86-P1/WT_Extra/20250805/WT/25805011.abf" 
),c(
  "base",
  "extra1.0",  "extra1.3",  "extra1.6",  "extra1.9",
  "extra2.0",  "extra2.3",  "extra2.6","extra2.9","extra2.12")),
nrow=2)

cellname<-c(paste("cell",c(2:5),sep=""))
cells<-list(cell2, cell3, cell4, cell5)
cell_values<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(cell_values)<-c("cell","condition","amplitude","v1/2","tau")
cell_values_act<-data.frame(matrix(ncol = 6, nrow = 0))
colnames(cell_values_act) <- c("cell", "condition", "voltage", "cond. norm.", "tail. curr.", "tau")


for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(rep(cellname[i],nrow(curr_cell)),
                       curr_cell[2,],
                       rep(0, nrow(curr_cell)),
                       rep(0, nrow(curr_cell)),
                       rep(0, nrow(curr_cell)))
  
  for (s in 1:ncol(curr_cell)){
    
    sweepnr<-11
    data<-readABF(curr_cell[1,s])
    cond<-rep(0, sweepnr)
    tail_curr<-rep(0, sweepnr)
    curr<-rep(0, sweepnr)
    volt<-10*(-6:4)
    tau <-rep(0, sweepnr)
    for (ii in 1:sweepnr){
      sweep.data<-as.data.frame(data,sweep=ii)
      sweep.max<-sweep.data[c(800:20000),c(1,4)]
      curr[ii]<-max(sweep.max)
      cond[ii]<-max(sweep.max/(volt[i]+95))
      sweep.tail<-sweep.data[c(45000:46000),c(1,4)]
      tail_curr[ii]<-min(sweep.tail)
      # --- Tau of Inactivation ---
      # Fit decay after peak
      peak_idx <- which(sweep.data[,4]==max(sweep.data[c(800:20000),c(1,4)]))
      decay_data <- sweep.data[peak_idx:45000,c(1,4)]
      #decay_data[,2]<-decay_data[,2]-min(decay_data[,2])
      colnames(decay_data)<-c("time","current")
      try({
        fit <- nlsLM(current ~ A * exp(time / tau) + C,
                     start = list(A = decay_data$current[1], tau = -0.5, C = min(decay_data$current)),
                     control = nls.lm.control(maxiter = 500), data=decay_data)
        tau[ii] <- coef(fit)["tau"]
        decay_data$fit<-predict(fit)
        if(ii==112){
          print(
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
                                   rep(curr_cell[2,s], sweepnr),
                                   volt,
                                   cond_norm,
                                   tail_norm,
                                   tau))
    #cell_values[s,4]<-coef(model)[2]
    cell_values_i[s,3]<-curr[11] #40 mV
    cell_values_i[s,4]<-tau[11] #40 mV
  }
  
  cell_values<-rbind(cell_values,cell_values_i)
}


cell_values[,3]<-as.numeric(cell_values[,3])

for ( i in cellname){
  
  cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra1")),4]<-cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra1")),3]/cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"base")),3]
  if (any((cell_values[,1]==i)&(str_detect(cell_values[,2],"extra1.12")))){
    cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra2")),4]<-cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra2")),3]/cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra1.12")),3]
  } else{
    cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra2")),4]<-cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra2")),3]/cell_values[(cell_values[,1]==i)&(str_detect(cell_values[,2],"extra1.9")),3]
  }
}

cell_values [c('Condition', 'Time')]<- str_split_fixed(cell_values$V2,"\\.",2)
cell_values<-cell_values[!cell_values$Condition=="base",]
colnames(cell_values)<-c("cell"  ,      "ident"   ,     "raw.Amp"   ,     "norm.Amp"   , "Tau",    "Condition", "Time")
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


ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2_Extra.svg", width = 3.5, height = 2)
cell_values$Time<-as.character(cell_values$Time)
cell_values$Time<-factor(cell_values$Time, levels = c("0","3","6","9","12"))

kcna2_tox<-lm(norm.Amp~Condition*Time, data = cell_values)
anova(kcna2_tox)
lsmeans(kcna2_tox, pairwise ~ Condition | Time, adjust = "tukey")

summary <- cell_values %>%
  group_by( Condition, Time) %>% 
  summarise(meanAmp = mean(norm.Amp),
            sdAmp = sd(norm.Amp),
            nAMp = n()) %>%
  mutate(SEMAmp = sdAmp/sqrt(nAMp))
cell_values$Tau[21]<-NA
#############
#Inaktivation
################
cell_values$Tau[cell_values$Tau<(-10)]<-NA

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

ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2_tau.svg", width = 3.5, height = 2)
kcna2_tau<-lm(Tau~Condition*Time, data = cell_values)
anova(kcna2_tau)
lsmeans(kcna2_tau, pairwise ~ Condition | Time, adjust = "tukey")

summary <- na.omit(cell_values) %>%
  group_by( Condition, Time) %>% 
  summarise(meanTau = mean(Tau),
            sdTau = sd(Tau),
            nTau = n()) %>%
  mutate(SEMTau = sdTau/sqrt(nTau))

##################
#Plot examples
###################

data1<-readABF(cell06[1,1])
data1<-as.data.frame(data1,sweep=12)
data2<-readABF(cell06[6,1])
data2<-as.data.frame(data2,sweep=12)
data3<-readABF(cell06[11,1])
data3<-as.data.frame(data3,sweep=12)
ggplot(data1,aes(`Time [s]`,`IN 0C [pA]`))+
  geom_line(color="#252525",size=2)+
  geom_line(data=data2,color="#3c5396",size=2)+
  geom_line(data=data3,color="#b5595f",size=2)+
  ylim(c(-4500,5000))+
  theme_prism(base_size = 14)+
  theme(axis.title.x=element_blank(),        axis.text.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),        axis.text.y=element_blank(),
        legend.position = "none")
ggsave(filename = "D:/Peter/Analysis/KCNA2/BMK86-P1/KCNA2_example.svg", width = 2, height = 2)
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
