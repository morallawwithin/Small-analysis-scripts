library(tidyverse)
library(cowplot)
#set your working directory
setwd("G:/Arbeit/EEG")
#now lets read all .txt files in the following folder (they should contain the exported raw values of the EEG)
filenames <- list.files("G:/Arbeit/EEG", pattern="*.txt", full.names=TRUE)
#this plots all files and saves a .svg image of them
lapply(filenames, function(x){
data<-read.table(file=x,sep = "",dec=",")
colnames(data)<-c("t","V")
data$time<-seq_along(data$t)
a<-ggplot(data,aes(x=time,y=V))+
  geom_line(lwd=1)+
  ylim(c(-4e-04,4e-04))+ #this limits the y axis from -400 to 400 µV
  theme_half_open()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_segment(aes(x = max(time)-max(time)/8, y = -3.5e-04, xend = max(time), yend = -3.5e-04),lwd=2)+ #the length of the x bar is 1/8 of the total length of the recording
  geom_segment(aes(x = max(time), y = -3.5e-04, xend = max(time), yend = -2.5e-04),lwd=2)              #the length of the y bar is 100 µV 
ggsave(plot=a,filename = gsub("txt$","svg",x))
})
