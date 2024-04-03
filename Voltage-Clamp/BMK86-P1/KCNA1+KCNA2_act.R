library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
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
                               pattern="_0(009[1-9]|101)"),
                           sep=""),
                    condition ), ncol = 2)
cell02<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(04[5-9]|5[0-5])"),
                           sep=""),
                     condition ), ncol = 2)
cell03<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(05[8-9]|6[0-8])"),
                           sep=""),
                     condition ), ncol = 2)
cell04<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(07[0-9]|80)"),
                           sep=""),
                     condition ), ncol = 2)
cell05<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_11",
                                      pattern="_00(08[2-9]|9[0-2])"),
                           sep=""),
                     condition ), ncol = 2)
cell06<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12",
                                      pattern="_00(00[1-9]|1[0-1])"),
                           sep=""),
                     condition ), ncol = 2)
cell07<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1+KCNA2/2022_4_12",
                                      pattern="_00(03[6-9]|4[0-6])"),
                           sep=""),
                     condition ), ncol = 2)

cell_values_act<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(cell_values_act)<-c("cell","voltage","cond. norm.","tail. curr.")
cell_values<-data.frame(matrix(ncol = 3, nrow = 0))
colnames(cell_values)<-c("cell","amplitude","v1/2")
cellname<-c(paste("cell0",c(1:7),sep=""))
cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07)
for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(cellname[i],0,0)
  s<-1
  sweepnr<-16
  data<-readABF(curr_cell[1,s])
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
  
  
  sweep<-data.frame("voltage"=volt, 
                    "conductance"= cond_norm,
                    "tail_curr"=tail_norm)
  SS<-getInitial(conductance~SSlogis(voltage,alpha,xmid,scale),data=sweep)
  
  activation <- function(g, Vhalf, k,c,V) (g/(1+exp((Vhalf-V)/k))+c)
  model <- nls(conductance ~ activation(myg,myVhalf,myk,myc,voltage), data=sweep, start=list(myg=SS["alpha"],myVhalf=SS["xmid"],myk=SS["scale"],myc=0),control = nls.control(maxiter = 400))
  
  curve_df <- data.frame("voltage" = seq(min(sweep$voltage), max(sweep$voltage), length.out = 100),
                         "conductance" = seq(min(sweep$conductance), max(sweep$conductance), length.out = 100))
  curve_df$prob <- predict(model, curve_df, type = "response")
  print(
    ggplot(curve_df, aes(x = voltage, y = prob)) +
      geom_line() +
      geom_point(data=sweep,aes(x=voltage,y=conductance))+
      geom_point(data=sweep,aes(x=voltage,y=tail_curr), col="blue")+
      theme_classic())
  
  cell_values_i[s,3]<-coef(model)[2]
  cell_values_i[s,2]<-curr[12]
  cell_values<-rbind(cell_values,cell_values_i)
  cell_values_act<-rbind(cell_values_act,
                         cbind(rep(cellname[i],16),
                               volt,
                               cond_norm,
                               tail_norm))
}
kcna1_kcna2<-as.numeric(cell_values[,3])
kcna1_kcna2_act<-cell_values_act
