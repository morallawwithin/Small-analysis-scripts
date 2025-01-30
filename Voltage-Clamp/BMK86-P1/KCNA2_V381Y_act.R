library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
cell1<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-241226/24d27000.abf"
cell2<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250106/25106002.abf"  
#cell3<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250106/25106007.abf" 
cell4<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120000.abf" 
cell5<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120005.abf" 
cell6<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120006.abf" 
cell7<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120007.abf" 
#cell8<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120008.abf"
cell9<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120010.abf"
cell10<-  "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0003.abf"
cell11<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0013a.abf"
#cell12<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0000.abf"
cell13<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0007.abf"
cell14<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0007.abf"
cell15<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0003_#1_Baseline.abf"
#cell16<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0014_#2_Baseline.abf"
cell17<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0012.abf"
cell_values<-data.frame(matrix(ncol = 3, nrow = 0))
colnames(cell_values)<-c("cell","amplitude","v1/2")
cell_values_act<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(cell_values_act)<-c("cell","voltage","cond. norm.","tail. curr.")
cells<-list(cell1,cell2,cell4,cell5,cell6,cell7,cell9,
            cell10,cell11,cell13,cell14,cell15,cell17)
cellname<-c("cell1","cell2","cell4","cell5","cell6","cell7","cell9",
            "cell10","cell11","cell13","cell14","cell15","cell17")
Ek<-(-75)

for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(cellname[i],0,0)
    data<-readABF(curr_cell)
    sweepnr<-length(data$data)
    cond<-rep(0, sweepnr)
    tail_curr<-rep(0, sweepnr)
    curr<-rep(0, sweepnr)
    GHK<-rep(0, sweepnr)
    if(sweepnr==16){
      volt<-10*(-8:7)+10^-20
      for (ii in 1:sweepnr){
      sweep.data<-as.data.frame(data,sweep=ii)
      sweep.max<-sweep.data[c(410:5000),c(1,4)]
      curr[ii]<-max(sweep.max)
      cond[ii]<-max(sweep.max/(volt[ii]-Ek))
      GHK[ii]<-max(sweep.max/
                     (volt[ii]/25*
                        (exp((volt[ii]-Ek)/25)-1)/
                        (exp((volt[ii])/25)-1)
                      ))
      sweep.tail<-sweep.data[c(20410:20500),c(1,4)]
      tail_curr[ii]<-min(sweep.tail)
    } 
    }else if(sweepnr==11){
      volt<-10*(-6:4)+10^-20
      for (ii in 1:sweepnr){
        sweep.data<-as.data.frame(data,sweep=ii)
        sweep.max<-sweep.data[c(500:5000),c(1,4)]
        curr[ii]<-max(sweep.max)
        cond[ii]<-max(sweep.max/(volt[ii]-Ek))
        GHK[ii]<-max(sweep.max/
                       (volt[ii]/25*
                          (exp((volt[ii]-Ek)/25)-1)/
                          (exp(volt[ii]/25)-1)
                       ))
        sweep.tail<-sweep.data[c(45780:46000),c(1,4)]
        tail_curr[ii]<-min(sweep.tail)
      }
      }
    
    
    
    cond_norm<-cond/max(cond)
    tail_norm<-tail_curr/min(tail_curr)
    GHK_norm<-GHK/max(GHK)
    
    
    sweep<-data.frame("voltage"=volt, 
                      "conductance"= cond_norm,
                      "tail_conductance"=tail_norm,
                      "GHK"=GHK_norm)
    if(i==2){
      sweep$tail_conductance[12:16]<-NA
    }else if(i==12){
      sweep$tail_conductance[9:11]<-NA
    }
    SS<-getInitial(tail_conductance~SSlogis(voltage,alpha,xmid,scale),data=sweep)
    
    activation <- function(g, Vhalf, k,c,V) (g/(1+exp((Vhalf-V)/k))+c)
    model <- nls(tail_conductance ~ activation(myg,myVhalf,myk,myc,voltage),
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
        geom_point(data=sweep,aes(x=voltage,y=GHK),col="blue")+
      theme_classic())
  cell_values_act<-rbind(cell_values_act,
                         cbind(rep(cellname[i],sweepnr),
                         volt,
                         cond_norm,
                         tail_norm))  
  cell_values_i[1,3]<-coef(model)[2]
  cell_values<-rbind(cell_values,cell_values_i)
  
}
V381Y<-as.numeric(cell_values[,3])
V381Y_act<-cell_values_act

