library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(pspline)
setwd("D:/Peter/Data/Scna1a/M-current")

#select the cell
cell01<-"G:/Peter/20240322 - Copia/24322009 - Copia - filtered - leak sub.abf"
cell02<-"G:/Peter/20240322 - Copia/24322011 - Copia - filtered - leak sub.abf"
cell03<-"G:/Peter/20240322 - Copia/24322014 - Copia - filtered - leak sub.abf"
cell04<-"G:/Peter/20240322 - Copia/24322016 - Copia - filtered - leak sub.abf"
cell05<-"G:/Peter/20240322 - Copia/24322018 - Copia - filtered - leak sub.abf"
cell06<-"G:/Peter/20240322 - Copia/24322021 - Copia - filtered - leak sub_2.abf"
cell07<-"G:/Peter/20240322 - Copia/24322022 - Copia - filtered - leak sub.abf"
cell08<-"G:/Peter/20240322 - Copia/24322024 - Copia - filtered - leak sub.abf"

cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08)
cellname<-c("cell1","cell2","cell3","cell4", "cell5", "cell6", "cell7", "cell8")

sweep<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(sweep)<-c("cell","voltage","current","holding+20")

for ( i in 1:length(cellname)){

#
curr_cell<-cells[[i]]
#How many sweeps
  sweepnr<-12
  samplerate<-100000 #/s
#read data  
  data<-readABF(curr_cell)
#prepare to store  
  meanMcurr<-rep(0, sweepnr)
  IV_Mcurr<-rep(0, sweepnr)
#which voltage at which step  
  volt<-10*(2:-9)
#do the analysis  
for (ii in 1:sweepnr){
    sweep.data<-as.data.frame(data,sweep=ii)
    sweep.tail<-sweep.data[c((8*samplerate):(8.25*samplerate)),]
    meanMcurr[ii]<-mean(sweep.tail$`IN 2 [pA]`)
    
    sweep.IV<-sweep.data[c((8.2*samplerate):(9.2*samplerate)),]
    sweep.IV$first_derivate<-predict(sm.spline(sweep.IV$`Time [s]`, sweep.IV$`IN 2 [pA]`), sweep.IV$`Time [s]`, 1)
    min_first_deriv<-min(sweep.IV$first_derivate[1:20000])
    min_first_deriv_index<-match(min_first_deriv,sweep.IV$first_derivate)
    max_first_deriv<-max(sweep.IV$first_derivate[min_first_deriv_index:20000])
    ten_percent<-min_first_deriv+0.9*(max_first_deriv-min_first_deriv)
    start_ind<-which(sweep.IV$first_derivate[min_first_deriv_index:20000]>ten_percent)[1]
    start<-start_ind+min_first_deriv_index
    start_value<-sweep.IV$`IN 2 [pA]`[start]
    end_value<-mean(sweep.IV$`IN 2 [pA]`[(0.7*samplerate):(0.85*samplerate)])
    
#    print(ggplot(sweep.IV,aes(`Time [s]`,`IN 2 [pA]`))+
#      geom_line() +
#      geom_point(aes(x=sweep.IV$`Time [s]`[start],y=start_value),col="red", size=2)+ #     
#      geom_line(aes(x=sweep.IV$`Time [s]`[c(0.7*samplerate,0.85*samplerate)],y=c(end_value,end_value)),col="red")+
#      theme_classic())
    
    IV_Mcurr[ii]<-start_value-end_value

}

  sweep<-rbind(sweep,
    data.frame("cell"=rep(cellname[i],sweepnr),
                    "voltage"=volt, 
                    "current"= IV_Mcurr,
                    "holding+20"=meanMcurr))

}    
ggplot(sweep,aes(voltage,current, col=cell))+
       geom_point()+
       theme_classic()
  
write_xlsx(sweep, path=paste0(curr_cell[1],"_Mcurrent.xlsx"))
  
  