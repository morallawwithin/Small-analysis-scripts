library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
cell1<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0003.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0004.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0005.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0006.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0007.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0008.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0009.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0010.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0011.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0012.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0013.abf","toxin.12"),
              nrow=2)

cell2<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0013a.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0014.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0015.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0016.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0017.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0018.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0019.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0020.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0021.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0022.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0023.abf","toxin.12"),
              nrow=2)

cell3<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0000.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0001.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0002.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0003.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0004.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0005.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0006.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0007.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0008.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0009.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0010.abf","toxin.12"),
              nrow=2)

cell4<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0007.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0008.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0009.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0010.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0011.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0012.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0013.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0014.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0015.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0016.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0017.abf","toxin.12"),
              nrow=2)

cell5<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0007.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0008.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0009.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0010.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0011.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0012.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0013.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0014.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0015.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0016.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0017.abf","toxin.12"),
              nrow=2)

cell6<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0003_#1_Baseline.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0004_#1_WOBMK86-p1.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0005_#1_WOBMK86-p1.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0006_#1_WOBMK86-p1.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0007_#1_WOBMK86-p1.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0008#1_1_BMK86-p1.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0009_#1_2_BMK86-p1.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0010_#1_3_BMK86-p1.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0011_#1_4_BMK86-p1.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0012_#1_5_BMK86-p1.abf","toxin.12"),
              nrow=2)

cell7<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0014_#2_Baseline.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0015_#2_WOBMK86-p1.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0016_#2_WOBMK86-p1.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0017_#2_WOBMK86-p1.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0018_#2_WOBMK86-p1.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0019_#2_BMK86-p1.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0020_#2_BMK86-p1.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0021_#2_BMK86-p1.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0022_#2_BMK86-p1.abf","toxin.9"),
              nrow=2)

cell8<-matrix(c("D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0012.abf","base",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0013.abf","extra.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0014.abf","extra.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0015.abf","extra.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0016.abf","extra.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0017.abf","extra.12",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0018.abf","toxin.0",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0019.abf","toxin.3",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0020.abf","toxin.6",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0021.abf","toxin.9",
                "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0022.abf","toxin.12"),
              nrow=2)
#curr_cell<-cell5
cell_values<-data.frame(matrix(ncol = 3, nrow = 0))
colnames(cell_values)<-c("cell","amplitude","v1/2")
cell_values_act<-data.frame(matrix(ncol = 3, nrow = 0))
colnames(cell_values)<-c("cell","voltage","cond. norm.")
cells<-list(cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8)
cellname<-c("cell1","cell2","cell3","cell4","cell5","cell6","cell7","cell8")
for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(cellname[i],0,0)
    s<-1
    data<-readABF(curr_cell[1,s])
    cond<-rep(0, 11)
    tail_curr<-rep(0, 11)
    curr<-rep(0, 11)
    volt<-(-6:4)*10
    
    for (ii in 1:11){
      sweep.data<-as.data.frame(data,sweep=ii)
      sweep.max<-sweep.data[c(780:45000),c(1,4)]
      curr[ii]<-max(sweep.max)
      cond[ii]<-max(sweep.max/(volt[i]+95))
      sweep.tail<-sweep.data[c(45000:46000),c(1,4)]
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
  cell_values_act<-rbind(cell_values_act,
                         cbind(rep(cellname[i],11),
                         volt,
                         cond_norm))  
  cell_values_i[s,3]<-coef(model)[2]
  cell_values<-rbind(cell_values,cell_values_i)
  
}
V381Y<-as.numeric(cell_values[,3])
V381Y_act<-cell_values_act
