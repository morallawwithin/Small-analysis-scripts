library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(pspline)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys")

####################
##select the cells
####################
directory<-"D:/Peter/Data/KCNA2/P405L_Mice/Ephys/Peter"
cell01<-c(paste0(directory,
                 "/20240503/24503002",
                 ".abf"),
          "Kcna2+/+")
cell02<-c(paste0(directory,
                 "/20240506/24506007",
                 ".abf"),
          "Kcna2+/+")
cell03<-c(paste0(directory,
                 "/20240506/24506015",
                 ".abf"),
          "Kcna2+/+")
cell04<-c(paste0(directory,
                 "/20240506/24506022",
                 ".abf"),
          "Kcna2+/+")
cell05<-c(paste0(directory,
                 "/20240506/24506029",
                 ".abf"),
          "Kcna2+/+")
cell06<-c(paste0(directory,
                 "/20240507/24507004",
                 ".abf"),
          "Kcna2+/+")
cell07<-c(paste0(directory,
                 "/20240510/24510002",
                 ".abf"),
          "Kcna2+/+")
cell08<-c(paste0(directory,
                 "/20240510/24510015",
                 ".abf"),
          "Kcna2+/+")
cell09<-c(paste0(directory,
                 "/20240510/24510021",
                 ".abf"),
          "Kcna2+/+")
cell10<-c(paste0(directory,
                 "/20240604/24604014",
                 ".abf"),
          "Kcna2+/+")
cell11<-c(paste0(directory,
                 "/20240604/24604030",
                 ".abf"),
          "Kcna2+/+")
cell12<-c(paste0(directory,
                 "/20240610/24610006",
                 ".abf"),
          "Kcna2+/P405L")
cell13<-c(paste0(directory,
                 "/20240610/24610013",
                 ".abf"),
          "Kcna2+/P405L")
cell14<-c(paste0(directory,
                 "/20240611/24611001",
                 ".abf"),
          "Kcna2+/P405L")
cell15<-c(paste0(directory,
                 "/20240611/24611009",
                 ".abf"),
          "Kcna2+/P405L")
cell16<-c(paste0(directory,
                 "/20240613/24613001",
                 ".abf"),
          "Kcna2+/P405L")
cell17<-c(paste0(directory,
                 "/20240613/24613014",
                 ".abf"),
          "Kcna2+/P405L")
cell18<-c(paste0(directory,
                 "/20240502/24502001",
                 ".abf"),
          "Kcna2+/P405L")
##########
#prepare everything for the loop
##########
cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, 
            cell11, cell12, cell13, cell14, cell15, cell16,cell17,cell18)
cellname<-c(paste0("cell0",c(1:9)),paste0("cell",c(10:18)))

sag_data<-data.frame(matrix(ncol = 7, nrow = 0))
colnames(sag_data)<-c("cell","genotype","current","sag")
pass_properties<-data.frame(matrix(ncol = 9, nrow = 0))
colnames(pass_properties)<-c("cell","genotype","resting_memb_pot","input_resis")
#Samplerate
samplerate<-100000 #/s

########
#loop starts here
#######
for ( i in 1:length(cellname)){
  #select the cell
  curr_cell<-cells[[i]]
  #read data of one cell  
  data<-readABF(curr_cell[1])
  #How many sweeps
  sweepnr<-length(data$data)
  curr<-(-10)*(1:sweepnr)
  sag<-rep(0, sweepnr)
  resting_memb_pot<-rep(0, sweepnr)
  input_resis<-rep(0, sweepnr)
  for (ii in 1:sweepnr){
    #read the sweep data
    sweep.data<-as.data.frame(data,sweep=ii)
    #get the resting potential
    resting_memb_pot[ii]<-mean(sweep.data$`INcc 0 [mV]`[1:(0.01*samplerate)])
    min_pot<-min(sweep.data$`INcc 0 [mV]`[1:(0.25*samplerate)])
    input_resis[ii]<-(resting_memb_pot[ii]-min_pot)/abs(curr[ii])*1000
    sag[ii]<-mean(sweep.data$`INcc 0 [mV]`[(0.45*samplerate):(0.5*samplerate)])-min_pot
  }
  sag_data<-rbind(sag_data,
               data.frame("cell"=rep(cellname[i],sweepnr),
                          "genotype"=rep(curr_cell[2],sweepnr),
                          "current"=curr, 
                          "sag"= sag))
  pass_properties<-rbind(pass_properties,
                  data.frame("cell"=cellname[i],
                             "genotype"=curr_cell[2],
                             "resting_memb_pot"=mean(resting_memb_pot), 
                             "input_resis"= mean(input_resis)))
}
p1<-ggplot(pass_properties,aes(genotype,resting_memb_pot, fill=genotype,col=genotype))+
  geom_boxplot()+
  geom_point()+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  theme_prism(base_size = 14)+
  ylim(c(-100,0))+
  xlab("genotype") + ylab("resting_memb_pot")
ggsave(p1,width = 4, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/resting_memb_pot.png")
p2<-ggplot(pass_properties,aes(genotype,input_resis, fill=genotype,col=genotype))+
  geom_boxplot()+
  geom_point()+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  theme_prism(base_size = 14)+
  ylim(c(0,400))+
  xlab("genotype") + ylab(expression(paste("input resistance [M",Omega,"]")))
ggsave(p2,width = 4, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/input_resis.png")
p3<-ggplot(sag_data,aes(current,sag,group=as.factor(genotype), col=as.factor(genotype),fill=as.factor(genotype)))+  
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 30,size=1,  position = position_dodge(width = 0.5)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'path',  size=1, position = position_dodge(width = 0.5), aes(col=as.factor(genotype))) +
  stat_summary(fun = mean,
               geom = 'point', size=5, position = position_dodge(width = 0.5),shape=17) +
  scale_colour_manual(values = c("black", "blue")) +
  ylim(c(0,3))+
  theme_prism(base_size = 14)+
  xlab("injected current [pA]") + ylab("sag potential [mV]")
ggsave(p3,width = 6, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/sag-pot.png")
