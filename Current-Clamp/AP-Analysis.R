library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(pspline)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys")

#select the cell
cell01<-"D:/Peter/Data/KCNA2/P405L_Mice/Ephys/Peter/20240610/24610015.abf"


cells<-list(cell01)#, cell02, cell03, cell04, cell05, cell06, cell07, cell08)
cellname<-c("cell1","cell2","cell3","cell4", "cell5", "cell6", "cell7", "cell8")

sweep<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(sweep)<-c("cell","current","AP")

AP_properties<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(AP_properties)<-c("cell","current","AP_Nr","Threshold")

AP_IFF<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(AP_IFF)<-c("cell","current","AP_Nr","IFF")

for ( i in 1:length(cellname)){

#
curr_cell<-cells[[i]]
#How many sweeps
  sweepnr<-17
  samplerate<-100000 #/s
#read data  
  data<-readABF(curr_cell)
#prepare to store  
  AP_nr<-rep(0, sweepnr)
  AP_thres<-list()
  AP_isi<-list()
#which voltage at which step  
  curr<-25*(-4:12)
#do the analysis  
for (ii in 1:sweepnr){
    sweep.data<-as.data.frame(data,sweep=ii)


    
    sweep.data$first_derivate<-predict(sm.spline(sweep.data$`Time [s]`, sweep.data$`INcc 0 [mV]`), sweep.data$`Time [s]`, 1)/samplerate*1000

    
    AP_ind_raw<-which(sweep.data$`INcc 0 [mV]`>0)#indices that are above 0 mV
    if(!is_empty(AP_ind_raw)){
      AP_ind_diff<-c(1000,diff(AP_ind_raw))#you loose one index this way, therefore you have to add them again
      Ap_ind_log<-AP_ind_diff>200 #larger difference than 2ms
      AP_ind_zero<-AP_ind_raw[Ap_ind_log]
      AP_ind<-c()
      for (iii in 1:length(AP_ind_zero)){
        AP_data<-sweep.data$first_derivate[(AP_ind_zero[iii]-500):(AP_ind_zero[iii])]#threshold last 3 ms
        AP_thres_ind<-which(AP_data>20)#indices that are above the firing threshold defined as 20 V/s
        AP_thres_diff<-c(2,diff(AP_thres_ind))
        AP_ind_thres<-AP_thres_ind[tail(which(AP_thres_diff>1),1)]
        
        AP_ind[iii]<-AP_ind_zero[iii]-500+AP_ind_thres #threshold
      }

      
      AP_nr[ii]<-length(AP_ind)
      AP_thres[[ii]]<-sweep.data$`INcc 0 [mV]`[AP_ind]    
      AP_isi[[ii]]<-1/diff(sweep.data$`Time [s]`[AP_ind] )
      
      
      AP_df<-data.frame('ind'=sweep.data$`Time [s]`[AP_ind],'volt'=sweep.data$`INcc 0 [mV]`[AP_ind])
      print(ggplot(sweep.data,aes(`Time [s]`,`INcc 0 [mV]`))+
      geom_line(aes(`Time [s]`,first_derivate), color='grey')+
      geom_line()+
      geom_point(data=AP_df,aes(`ind`,volt),color="red")+
      ylim(c(-80,50)))
    }
    else {
      AP_nr[ii]<-0
    }#+

    


}

  sweep<-rbind(sweep,
    data.frame("cell"=rep(cellname[i],sweepnr),
               "current"=curr, 
                    "AP"= AP_nr))
  
  
  names(AP_thres)<-curr
  AP_thres_names<-as.character( names(unlist(AP_thres)))
    AP_thres_nr<-sub("[123]?[2570][05]","",AP_thres_names)
  AP_thres_nr[AP_thres_nr=='']<-1
  
  AP_properties<-rbind(AP_properties,
                       data.frame(
                         "cell"=rep(cellname[i],length(unlist(AP_thres))),
                         "current"=as.numeric( stringr::str_extract(AP_thres_names,"[123]?[2570][05]")),
                         "AP_Nr"=as.numeric(AP_thres_nr),
                         "Threshold"=unlist(AP_thres)
                         ))
  names(AP_isi)<-curr
  AP_isi_names<-as.character( names(unlist(AP_isi)))
  AP_isi_nr<-sub("[123]?[2570][05]","",AP_isi_names)
  AP_isi_nr[AP_isi_nr=='']<-1
  
  AP_IFF<-rbind(AP_IFF,
                       data.frame(
                         "cell"=rep(cellname[i],length(unlist(AP_isi))),
                         "current"=as.numeric( stringr::str_extract(AP_isi_names,"[123]?[2570][05]")),
                         "AP_Nr"=as.numeric(AP_isi_nr),
                         "IFF"=unlist(AP_isi)
                       ))
  
}    
ggplot(AP_IFF,aes(AP_Nr,IFF, col=current))+
       geom_point()+
       theme_classic()

ggplot(AP_properties,aes(AP_Nr,Threshold, col=current))+
  geom_point()+
  theme_classic()
  
write_xlsx(sweep, path=paste0(curr_cell[1],"_Mcurrent.xlsx"))
  
  