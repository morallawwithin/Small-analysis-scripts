library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(pspline)
library(readxl)
library(ggbeeswarm)
library(scattermore)
library(lsmeans)
library(lme4)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys")

####################
##select the cells
####################
dataset<-"EC_L5PN"#"Cortex_L2&3_PN"
data <- read_excel(paste0(dataset,".xlsx"))
setwd(paste0("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/",dataset))
data<-data[data$protocol=="AP_Input",]
cells<-data[,c(2,3)]
cellname<-data$cell
##########
#prepare everything for the loop
##########
sag_data<-data.frame(matrix(ncol = 9, nrow = 0))
colnames(sag_data)<-c("cell","genotype","current","sag","peak_sag","steady_state")
pass_properties<-data.frame(matrix(ncol = 9, nrow = 0))
colnames(pass_properties)<-c("cell","genotype","resting_memb_pot","input_resis")
#Samplerate
samplerate<-100000 #/s

########
#loop starts here
#######
for ( i in 1:length(cellname)){
  #select the cell
  curr_cell<-cells[i,]
  #read data of one cell  
  data<-readABF(curr_cell$file)
  #How many sweeps
  sweepnr<-length(data$data)
  curr<-(-10)*(1:sweepnr)
  sag<-rep(0, sweepnr)
  resting_memb_pot<-rep(0, sweepnr)
  input_resis<-rep(0, sweepnr)
  peak_sag<-rep(0, sweepnr)
  steady_state<-rep(0, sweepnr)
  for (ii in 1:sweepnr){
    #read the sweep data
    sweep.data<-as.data.frame(data,sweep=ii)
    #get the resting potential
    resting_memb_pot[ii]<-mean(sweep.data[,2][1:(0.01*samplerate)])
    min_pot<-min(sweep.data[,2][1:(0.25*samplerate)])
    input_resis[ii]<-(resting_memb_pot[ii]-min_pot)/abs(curr[ii])*1000
    sag[ii]<-mean(sweep.data[,2][(0.45*samplerate):(0.5*samplerate)])-min_pot
    steady_state[ii]<-mean(sweep.data[,2][(0.45*samplerate):(0.5*samplerate)])
    peak_sag[ii]<-min_pot
  }
  sag_data<-rbind(sag_data,
               data.frame("cell"=rep(cellname[i],sweepnr),
                          "genotype"=rep(curr_cell$genotype,sweepnr),
                          "current"=curr, 
                          "sag"= sag,
                          "peak_sag"=peak_sag,
                          "steady_state"=steady_state))
  pass_properties<-rbind(pass_properties,
                  data.frame("cell"=cellname[i],
                             "genotype"=curr_cell$genotype,
                             "resting_memb_pot"=mean(resting_memb_pot), 
                             "input_resis"= mean(input_resis)))
}
p1<-ggplot(pass_properties,aes(genotype,resting_memb_pot, fill=genotype,col=genotype))+
  geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  theme_prism(base_size = 14)+
  ylim(c(-100,-50))+
  xlab("genotype") + ylab("rest. memb. pot.")+
  theme(legend.position = "none") 
ggsave(p1,width = 3, height = 4,
       file="resting_memb.png")
p2<-ggplot(pass_properties,aes(genotype,input_resis, fill=genotype,col=genotype))+
  geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  theme_prism(base_size = 14)+
  ylim(c(0,600))+
  xlab("genotype") + ylab(expression(paste("input resistance [M",Omega,"]")))+
  theme(legend.position = "none")
ggsave(p2,width = 3, height = 4,
       file="input_resis.png")
p3<-ggplot(sag_data,aes(current,sag,group=as.factor(genotype), col=as.factor(genotype),fill=as.factor(genotype)))+  
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 10,size=1,  position = position_dodge(width = 0.5)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'path',  size=1, position = position_dodge(width = 0.5), aes(col=as.factor(genotype))) +
  stat_summary(fun = mean,
               geom = 'point', size=5, position = position_dodge(width = 0.5),shape=17) +
  scale_colour_manual(values = c("black", "blue")) +
  #ylim(c(-0.5,6))+
  theme_prism(base_size = 14)+
  xlab("injected current [pA]") + ylab("sag potential [mV]")
ggsave(p3,width = 6, height = 4,
       file="sag-pot.png")

p5<-ggplot(sag_data,aes(steady_state,peak_sag, col=as.factor(genotype),fill=as.factor(genotype)))+
  geom_point(shape=21, size=4)+
  geom_smooth(method='lm', formula= y~x)+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  theme_prism(base_size = 14)+
  ylim(c(-50,-120))+
  xlim(c(-50,-120))+  
  xlab("steady state potential [mV]") + ylab("peak sag potential [mV]")+
  theme(legend.position = "none")
p5
ggsave(p5,width = 6, height = 4,
       file="sag-pot_vs_memb.png")
sag.mdl<-lm(peak_sag~steady_state*genotype,data=sag_data)
summary(sag.mdl)

data1<-readABF(cells$file[17])#17#1
data12<-as.data.frame(data1,sweep=7)

data2<-readABF(cells$file[9])#9#22
data22<-as.data.frame(data2,sweep=9)


p4<-ggplot(data12,aes(data12[,1],data12[,2]))+
  geom_line()+
  geom_line(data=data22,aes(data22[,1],data22[,2]),color="blue")+
  #geom_scattermore(shape='.')+
  #geom_scattermore(data=data22,aes(data22[,1],data22[,2]),color="blue",shape='.')+
  ylim(c(-105,-65))+
  theme_prism(base_size = 14)+
  xlab("time [s]") + ylab("memb. pot. [mV]")+
  geom_segment(aes(x = 0.7, y = -95, xend = 0.9, yend = -95),lwd=2,col= "black")+ #the length of the x bar is 1/8 of the total length of the recording
  geom_segment(aes(x = 0.9, y = -95, xend = 0.9, yend = -85),lwd=2,col= "black")   +           #the length of the y bar is 100 µV 
  theme_prism(base_size = 14)+
  geom_text(x=0.96, y=-90, label="10 mV",col= "black")+
  geom_text(x=0.8, y=-96.5, label="200 ms",col= "black")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
p4
ggsave(p4,width = 6,height = 4,
       file="examples_sag.png")  
