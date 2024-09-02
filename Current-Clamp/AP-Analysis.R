library(readABF)
library(tidyverse)
library(readxl)
library(ggprism)
library(pspline)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys")

####################
##select the cells
####################
dataset<-"Cortex_L2&3_PN"
data <- read_excel(paste0(dataset,".xlsx"))
setwd(paste0("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/",dataset))
data<-data[data$protocol=="AP",]
cells<-data[,c(2,3)]
cellname<-data$cell

##########
#prepare everything for the loop
##########
sweep<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(sweep)<-c("cell","genotype","current","AP")
AP_properties<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(AP_properties)<-c("cell","genotype","current","AP_Nr","Threshold")
AP_IFF<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(AP_IFF)<-c("cell","genotype","current","AP_Nr","IFF")

########
#loop starts here
#######
for ( i in 1:length(cellname)){
#select the cell
  curr_cell<-cells[i,]
#read data of one cell  
  data<-readABF(curr_cell$file)
#How many sweeps
  sweepnr<-length(data$data)
  samplerate<-100000 #/s
#prepare to store AP parameters temporarly 
  AP_nr<-rep(0, sweepnr)
  AP_thres<-list()
  AP_isi<-list()
#which voltage at which step  
  curr<-25*(-4:(sweepnr-5))  
    #do the sweep analysis  
    for (ii in 1:sweepnr){
      #read the sweep data
      sweep.data<-as.data.frame(data,sweep=ii) 
      #calculate the first derivative
      sweep.data$first_derivate<-predict(sm.spline(sweep.data[,1], sweep.data[,2]), sweep.data[,1], 1)/samplerate*1000
      #Find action potentials
      #first: find indices that are above 0 mV
      AP_ind_raw<-which(sweep.data[,2]>0)
      #if there are any continue analysis
      if(!is_empty(AP_ind_raw)){
        #calculate the difference between the indices >0 mV
        #you loose one index this way, therefore you have to add them again
        AP_ind_diff<-c(1000,diff(AP_ind_raw))
        #select difference larger than 2ms
        Ap_ind_log<-AP_ind_diff>200 
        AP_ind_zero<-AP_ind_raw[Ap_ind_log]
        #prepare for 2. criterion which also finds thresholds
        AP_ind<-c()
        for (iii in 1:length(AP_ind_zero)){
          #take data 5ms before first value above 0mV
          AP_data<-sweep.data$first_derivate[(AP_ind_zero[iii]-500):(AP_ind_zero[iii])]
          #select indices that are above the firing threshold defined as 20 V/s
          AP_thres_ind<-which(AP_data>20)
          #calculate the difference between the indices >20 V/s and add one to the start
          AP_thres_diff<-c(2,diff(AP_thres_ind))
          #select the last index of the ones where thedifference to the previous index was >1
          AP_ind_thres<-AP_thres_ind[tail(which(AP_thres_diff>1),1)]
          #store this index
          AP_ind[iii]<-AP_ind_zero[iii]-500+AP_ind_thres #threshold
      }
        #store the number of APs
        AP_nr[ii]<-length(AP_ind)
        #store the threshold of APs
        AP_thres[[ii]]<-sweep.data[,2][AP_ind]
        #store the IFF of APs
        AP_isi[[ii]]<-1/diff(sweep.data[,1][AP_ind] )
###############      
#uncomment the section below to display thresolds for each analyzed sweep
#############        
      #AP_df<-data.frame('ind'=sweep.data$`Time [s]`[AP_ind],'volt'=sweep.data$`INcc 0 [mV]`[AP_ind])
      #print(
      #ggplot(sweep.data,aes(`Time [s]`,`INcc 0 [mV]`))+
      #geom_line(aes(`Time [s]`,first_derivate), color='grey')+
      #geom_line()+
      #geom_point(data=AP_df,aes(`ind`,volt),color="red")+
      #ylim(c(-80,50))
      #)
    }
    else {
      AP_nr[ii]<-0
    }}

  #store the number of APs per cell
  sweep<-rbind(sweep,
    data.frame("cell"=rep(cellname[i],sweepnr),
               "genotype"=rep(curr_cell$genotype,sweepnr),
               "current"=curr[1:sweepnr], 
                    "AP"= AP_nr))
  
  #store the threshold of APs per cell and associate the number of the spike to the threshold
  names(AP_thres)<-curr
  AP_thres_names<-as.character( names(unlist(AP_thres)))
  AP_thres_nr<-sub("[123]?[2570][05]","",AP_thres_names)
  AP_thres_nr[AP_thres_nr=='']<-1
  AP_properties<-rbind(AP_properties,
                       data.frame(
                         "cell"=rep(cellname[i],length(unlist(AP_thres))),
                         "genotype"=rep(curr_cell$genotype,length(unlist(AP_thres))),
                         "current"=as.numeric( stringr::str_extract(AP_thres_names,"[123]?[2570][05]")),
                         "AP_Nr"=as.numeric(AP_thres_nr),
                         "Threshold"=unlist(AP_thres)
                         ))
  #store the threshold of APs per cell and associate the number of the spike to the threshold
  names(AP_isi)<-curr
  AP_isi_names<-as.character( names(unlist(AP_isi)))
  AP_isi_nr<-sub("[123]?[2570][05]","",AP_isi_names)
  AP_isi_nr[AP_isi_nr=='']<-1
  AP_IFF<-rbind(AP_IFF,
                       data.frame(
                         "cell"=rep(cellname[i],length(unlist(AP_isi))),
                         "genotype"=rep(curr_cell$genotype,length(unlist(AP_isi))),
                         "current"=as.numeric( stringr::str_extract(AP_isi_names,"[123]?[2570][05]")),
                         "AP_Nr"=as.numeric(AP_isi_nr),
                         "IFF"=unlist(AP_isi)
                       ))
 print(paste0("finished analysis of ",curr_cell$file)) 
}    

p1<-ggplot(sweep,aes(current,AP,group=as.factor(genotype), col=as.factor(genotype),fill=as.factor(genotype)))+  
  stat_summary(fun = mean, 
               fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = 'errorbar',  width = 30,size=1,  position = position_dodge(width = 0.5)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = 'path',  size=1, position = position_dodge(width = 0.5), aes(col=as.factor(genotype))) +
  stat_summary(fun = mean,
               geom = 'point', size=5, position = position_dodge(width = 0.5),shape=17) +
  scale_colour_manual(values = c("black", "blue")) +
  theme_prism(base_size = 14)+
  xlab("injected current [pA]") + ylab("number of AP")

model_AP<-lm(AP~current*genotype,data= sweep)
summary(model_AP)
ggsave(p1,width = 6, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/APnumber.png")

p2<-ggplot(AP_IFF,aes(as.factor(AP_Nr),IFF, col=genotype))+
  geom_boxplot()+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1))) +
  theme_prism(base_size = 14)+
  xlab("AP number") + ylab("inst. firing fre. [Hz]")
ggsave(p2,width = 6, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/instfiringfreqbox.png")

model_IFF<-lm(IFF~AP_Nr+current+genotype,data= AP_IFF)
summary(model_IFF)


ggplot(AP_properties,aes(as.factor(AP_Nr),Threshold, col=genotype))+
  geom_boxplot()+
  #scale_colour_manual(values = c("black", "blue")) +
  theme_prism(base_size = 14)+
  xlab("AP number") + ylab("Threshold")
  

##Examples
data1<-readABF(cell01[1])
data1<-as.data.frame(data1,sweep=13)
data2<-readABF(cell14[1])
data2<-as.data.frame(data2,sweep=13)
p3<-ggplot(data1,aes(`Time [s]`,`INcc 0 [mV]`))+
  geom_line()+
  geom_line(data=data2,color="blue")+
  ylim(c(-80,50))+
  theme_prism(base_size = 14)+
  xlab("time [s]") + ylab("memb. pot. [mV]")
ggsave(p3,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/examples_200pA.png")  
