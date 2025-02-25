library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
#cell1<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-241226/24d27000.abf" # 15%leak
cell2<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250106/25106002.abf"  
#cell3<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250106/25106007.abf" 
cell4<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120000.abf" 
cell5<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120005.abf" #?
cell6<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120006.abf" 
#cell7<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120007.abf" 
#cell8<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120008.abf"
cell9<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA_V381/KCNA2-250120/25120010.abf"
 cell10<-  "D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0003.abf"
 cell11<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_01_19_0013a.abf"#slow
 #cell12<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_02_0000.abf" #slow
 cell13<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_02_28_0007.abf"
 cell14<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_03_10_0007.abf"
 cell15<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0003_#1_Baseline.abf"
 #cell16<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_05_30_0014_#2_Baseline.abf"
 #cell17<-"D:/Peter/Data/KCNA2/BMK86/CHO_KCNA2_V381Y_BMK86-P1/2023_06_20_0012.abf"

cell_values<-data.frame(matrix(ncol = 3, nrow = 0))
colnames(cell_values)<-c("cell","amplitude","v1/2")
cell_values_act<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(cell_values_act)<-c("cell","voltage","cond. norm.","tail. curr.")
cells<-list(cell2,cell4,cell5,cell6,cell9,
            cell10,cell11,cell13,cell14,cell15,cell17)
cellname<-c("cell2","cell4","cell5","cell6","cell9",
            "cell10","cell11","cell13","cell14","cell15")
Ek<-(-75)

for ( i in 1:length(cellname)){
  curr_cell<-cells[[i]]
  cell_values_i<-cbind(cellname[i],0,0,0,0,0)
  data<-readABF(curr_cell)
  sweepnr<-length(data$data)
  curr<-rep(0, sweepnr)
  tail_curr<-rep(0, sweepnr)
  cond<-rep(0, sweepnr)
    if(sweepnr==16){
      volt<-10*(-8:7)
      for (ii in 2:sweepnr){
      sweep.data<-as.data.frame(data,sweep=ii)
      sweep.max<-sweep.data[c(420:20000),c(1,2)]
      curr[ii]<-max(sweep.max)
      cond[ii]<-max(sweep.max/(volt[ii]-Ek))
      sweep.tail<-sweep.data[c(20410:20500),c(1,4)]
      tail_curr[ii]<-min(sweep.tail)
      } 
      volt<-10*(-7:7)
      tail_curr<-tail_curr[2:sweepnr]
      cond<-cond[2:sweepnr]
    }else if(sweepnr==11){
      volt<-10*(-6:4)
      for (ii in 1:sweepnr){
        sweep.data<-as.data.frame(data,sweep=ii)
        sweep.max<-sweep.data[c(800:20000),c(1,2)]
        curr[ii]<-max(sweep.max)
        cond[ii]<-max(sweep.max/(volt[ii]-Ek))
          sweep.tail<-sweep.data[c(45780:46000),c(1,4)]
        tail_curr[ii]<-min(sweep.tail)
      }
      }
    
    
    
    cond_norm<-cond/max(cond)
    tail_norm<-tail_curr/min(tail_curr)

    
    sweep<-data.frame("voltage"=volt, 
                      "conductance"= cond_norm,
                      "tail_conductance"=tail_norm)
    #sweep<-sweep[1:which(tail_norm==1),]
    #SS<-getInitial(tail_conductance~SSlogis(voltage,alpha,xmid,scale),data=sweep)
    
    activation <- function(g, Vhalf, k,c,V) (g/(1+exp((Vhalf-V)/k))+c)
    model <- nls(conductance ~ activation(myg,myVhalf,myk,myc,voltage),
                 data=sweep, 
                 #start=list(myg=SS["alpha"],myVhalf=SS["xmid"],myk=SS["scale"],myc=0),
                start=list(myg=1,myVhalf=-5,myk=9,myc=0),
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
  cell_values_act<-rbind(cell_values_act,
                         cbind(rep(cellname[i],length(sweep$voltage)),
                               sweep))  
  cell_values_i[1,3:6]<-coef(model)
  cell_values_i[1,2]<-curr[12]
  cell_values<-rbind(cell_values,cell_values_i)
  
}

V381Y<-cell_values
V381Y_act<-cell_values_act
colnames(V381Y_act)<-c("cell","voltage","cond_norm.","tail_norm")
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