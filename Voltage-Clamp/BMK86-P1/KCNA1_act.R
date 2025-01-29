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

cell01<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                               pattern="_00(0[0-9]|10)"),
                           sep=""),
                    condition ), ncol = 2)
cell02<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_00(1[7-9]|2[0-7])"),
                           sep=""),
                     condition ), ncol = 2)
cell03<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_00(29|3[0-9])"),
                           sep=""),
                     condition ), ncol = 2)
cell04<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_00(4[0-9]|50)"),
                           sep=""),
                     condition ), ncol = 2)
cell05<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_00(5[2-9]|6[0-2])"),
                           sep=""),
                     condition ), ncol = 2)
cell06<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_00(6[3-9]|7[0-3])"),
                           sep=""),
                     condition ), ncol = 2)
cell07<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_00(7[5-9]|8[0-5])"),
                           sep=""),
                     condition ), ncol = 2)
cell08<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_01(1[2-9]|2[0-2])"),
                           sep=""),
                     condition ), ncol = 2)
cell09<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_7",
                                      pattern="_01(2[4-9]|3[0-4])"),
                           sep=""),
                     condition ), ncol = 2)
cell10<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_8/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_8",
                                      pattern="_00(0[0-9]|10)"),
                           sep=""),
                     condition ), ncol = 2)
cell11<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_8/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA1 WT/2022_4_8",
                                      pattern="_00(2[2-9]|3[0-2])"),
                           sep=""),
                     condition ), ncol = 2)
cell_values_act<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(cell_values_act)<-c("cell","voltage","cond. norm.","tail. curr.")
cell_values<-data.frame(matrix(ncol = 3, nrow = 0))
colnames(cell_values)<-c("cell","amplitude","v1/2")
cellname<-c(paste("cell0",c(1:9),sep=""),paste("cell",c(10:11),sep=""))
cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, cell11)
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
    cond[ii]<-max(sweep.max/(volt[ii]+95))
    sweep.tail<-sweep.data[c(20000:21000),c(1,3)]
    tail_curr[ii]<-min(sweep.tail)
  }
  cond_norm<-cond/max(cond)
  tail_norm<-tail_curr/min(tail_curr)
  
  
  sweep<-data.frame("voltage"=volt, 
                    "conductance"= cond_norm)
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
      theme_classic())
  
  cell_values_i[s,3]<-coef(model)[2]
  cell_values_i[s,2]<-curr[12]
  cell_values<-rbind(cell_values,cell_values_i)
  cell_values_act<-rbind(cell_values_act,
                         cbind(rep(cellname[i],11),
                               volt,
                               cond_norm,
                               tail_norm))
}
kcna1<-as.numeric(cell_values[,3])
kcna1_act<-cell_values_act
