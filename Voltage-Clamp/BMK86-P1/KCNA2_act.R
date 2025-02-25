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
# cell14<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11/",
#                            list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11",
#                                pattern="_00(2[0-9]|30)"),
#                            sep=""),
#                     condition ), ncol = 2)
cell14<-matrix(cbind(paste("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11/",
                           list.files("D:/Peter/Data/KCNA2/BMK86/KCNA2 WT/2022_4_11",
                               pattern="_01(1[4-9]|2[0-4])"),
                           sep=""),
                    condition ), ncol = 2)

cell_values<-data.frame(matrix(ncol = 3, nrow = 0))
colnames(cell_values)<-c("cell","amplitude","v1/2")
cell_values_act<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(cell_values_act)<-c("cell","voltage","cond. norm.","tail. curr.")
cellname<-c(paste("cell0",c(1:9),sep=""),paste("cell",c(10:14),sep=""))
cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, cell11, cell12, cell13, cell14)
for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(cellname[i],0,0,0,0,0)
  s<-1
  sweepnr<-16
  data<-readABF(curr_cell[1,s])
  cond<-rep(0, sweepnr)
  tail_curr<-rep(0, sweepnr)
  curr<-rep(0, sweepnr)
  volt<-10*(-8:7)
  
  for (ii in 1:sweepnr){
    sweep.data<-as.data.frame(data,sweep=ii)
    sweep.max<-sweep.data[c(500:10000),c(1,2)]
    curr[ii]<-max(sweep.max)
    cond[ii]<-max(sweep.max/(volt[ii]+75))
    sweep.tail<-sweep.data[c(20000:21000),c(1,3)]
    tail_curr[ii]<-min(sweep.tail)
  }
  volt<-10*(-7:7)
  tail_curr<-tail_curr[2:sweepnr]
  cond<-cond[2:sweepnr]
  cond_norm<-cond/max(cond)
  if(i==14){
  tail_curr[14:15]<-0  
  }
  
  tail_norm<-tail_curr/min(tail_curr)
  
  
  sweep<-data.frame("voltage"=volt, 
                    "conductance"= cond_norm,
                    "tail_conductance"=tail_norm)
  #sweep<-sweep[1:which(tail_norm==1),]
  SS<-getInitial(conductance~SSlogis(voltage,alpha,xmid,scale),data=sweep)
  
  activation <- function(g, Vhalf, k,c,V) (g/(1+exp((Vhalf-V)/k))+c)
  model <- nls(conductance ~ activation(myg,myVhalf,myk,myc,voltage),
               data=sweep, 
               start=list(myg=SS["alpha"],myVhalf=SS["xmid"],myk=SS["scale"],myc=0),
               #start=list(myg=1,myVhalf=0,myk=10,myc=0),
               control = nls.control(maxiter = 400))
  
  curve_df <- data.frame("voltage" = seq(min(sweep$voltage), max(sweep$voltage), length.out = 100),
                         "tail_conductance" = seq(min(tail_norm), max(tail_norm), length.out = 100))
  curve_df$prob <- predict(model, curve_df, type = "response")
  print(
    ggplot(curve_df, aes(x = voltage, y = prob)) +
      geom_line() +
      geom_point(data=sweep,aes(x=voltage,y=tail_conductance))+
      geom_point(data=sweep,aes(x=voltage,y=conductance),col="red")+
      theme_classic())
  
  cell_values_i[s,3:6]<-coef(model)
  cell_values_i[s,2]<-curr[12]
  cell_values<-rbind(cell_values,cell_values_i)
  cell_values_act<-rbind(cell_values_act,
                         cbind(rep(cellname[i],length(sweep$voltage)),
                               sweep))
}
wt<-cell_values
wt_act<-cell_values_act
colnames(wt_act)<-c("cell","voltage","cond_norm.","tail_norm")
