library(readABF)
library(tidyverse)
library(writexl)
library(ggprism)
library(pspline)
library(ggsignif)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys")

####################
##select the cells
####################
directory<-"D:/Peter/Data/KCNA2/P405L_Mice/Ephys/Peter"
directory2<-"D:/Peter/Data/KCNA2/P405L_Mice/Ephys/P405L_Liz/Ctx_SSP"
cell01<-c(paste0(directory,
                 "/20240503/24503003",
                 ".abf"),
          "Kcna2+/+")
cell02<-c(paste0(directory,
                 "/20240506/24506008",
                 ".abf"),
          "Kcna2+/+")
cell03<-c(paste0(directory,
                 "/20240506/24506016",
                 ".abf"),
          "Kcna2+/+")
cell04<-c(paste0(directory,
                 "/20240506/24506023",
                 ".abf"),
          "Kcna2+/+")
cell05<-c(paste0(directory,
                 "/20240506/24506030",
                 ".abf"),
          "Kcna2+/+")
cell06<-c(paste0(directory,
                 "/20240507/24507005",
                 ".abf"),
          "Kcna2+/+")
cell07<-c(paste0(directory,
                 "/20240510/24510003",
                 ".abf"),
          "Kcna2+/+")
cell08<-c(paste0(directory,
                 "/20240510/24510016",
                 ".abf"),
          "Kcna2+/+")
# cell09<-c(paste0(directory,
#                  "/20240510/24510023",
#                  ".abf"),
#           "Kcna2+/+")
cell10<-c(paste0(directory,
                 "/20240604/24604016",
                 ".abf"),
          "Kcna2+/+")
cell11<-c(paste0(directory,
                 "/20240604/24604032",
                 ".abf"),
          "Kcna2+/+")
cell12<-c(paste0(directory,
                 "/20240610/24610007",
                 ".abf"),
          "Kcna2+/P405L")
cell13<-c(paste0(directory,
                 "/20240610/24610014",
                 ".abf"),
          "Kcna2+/P405L")
cell14<-c(paste0(directory,
                 "/20240611/24611003",
                 ".abf"),
          "Kcna2+/P405L")
cell15<-c(paste0(directory,
                 "/20240611/24611010",
                 ".abf"),
          "Kcna2+/P405L")
cell16<-c(paste0(directory,
                 "/20240613/24613002",
                 ".abf"),
          "Kcna2+/P405L")
cell17<-c(paste0(directory,
                 "/20240613/24613015",
                 ".abf"),
          "Kcna2+/P405L")
cell18<-c(paste0(directory,
                 "/20240502/24502003",
                 ".abf"),
          "Kcna2+/P405L")
# cell19<-c(paste0(directory2,
#                  "/20240613/24613007",
#                  ".abf"),
#           "Kcna2+/P405L_Liz")
# cell20<-c(paste0(directory2,
#                  "/20240613/24613012",
#                  ".abf"),
#           "Kcna2+/P405L_Liz")
# cell21<-c(paste0(directory2,
#                  "/20240613/24613026",
#                  ".abf"),
#           "Kcna2+/P405L_Liz")
# cell22<-c(paste0(directory2,
#                  "/20240613/24613032",
#                  ".abf"),
#           "Kcna2+/P405L_Liz")
# cell23<-c(paste0(directory2,
#                  "/20240614/24614025",
#                  ".abf"),
#           "Kcna2+/P405L_Liz")
# cell24<-c(paste0(directory2,
#                  "/20240614/24614030",
#                  ".abf"),
#           "Kcna2+/P405L_Liz")
# cell25<-c(paste0(directory2,
#                  "/20240614/24614035",
#                  ".abf"),
#           "Kcna2+/P405L_Liz")
##########
#prepare everything for the loop
##########
cells<-list(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, #cell09, 
            cell10,            cell11, cell12, cell13, cell14, cell15, cell16,cell17,cell18#,cell19, cell20,
            #cell21, cell22, cell23, cell24, cell25
            )
cellname<-c("cell01" ,"cell02" ,"cell03" ,"cell04" ,"cell05" ,"cell06" ,"cell07" ,"cell08" ,#"cell09" ,
            paste0("cell",c(10:18))
            )
  #c(paste0("cell0",c(1:9)),paste0("cell",c(10:18)))

AP_data<-data.frame(matrix(ncol = 7, nrow = 0))
colnames(AP_data)<-c("volt","curr","time","ind","first_deriv","cell","genotype")
AP_properties<-data.frame(matrix(ncol = 9, nrow = 0))
colnames(AP_properties)<-c("cell","genotype","rheobase","threshold","FWHA","fAHP","risingetime","repolarizingtime","amplitude")
#Samplerate
samplerate<-100000 #/s

########
#loop starts here
#######
for ( i in 1:length(cellname)){

  #select the cell
  curr_cell<-cells[[i]]
  #read data of one cell  
  data<-readABF(curr_cell[1])
  sweep_length<-length(data$data[[1]][,1])/samplerate
  if (sweep_length<1.4){
    command_curr  <-c(rep(0,1875),c(1:samplerate)/(samplerate/200),rep(0,(0.2*samplerate-1875)))
  }else{
    command_curr  <-c(rep(0,2656),c(1:(1.5*samplerate))/((1.5*samplerate)/300),rep(0,(0.2*samplerate-2656)))
  }
  #plot(data$data[[1]][,2]*1000,type="l")
  #lines(command_curr,col="red")
  
  #find first AP
  AP_ind_raw<-head(which(data$data[[1]][,1]>0),1)
  range_after<-3500
  AP_data_cell<- data.frame(
    "volt" = c(data$data[[1]][(AP_ind_raw-600):(AP_ind_raw+range_after),1]),
    "curr" = command_curr[(AP_ind_raw-600):(AP_ind_raw+range_after)],
    "time" = (c(1:(sweep_length*samplerate))/samplerate)[(AP_ind_raw-600):(AP_ind_raw+range_after)],
    "ind" = c(1:(601+range_after)))
  AP_data_cell$first_deriv<-predict(sm.spline(AP_data_cell$time, AP_data_cell$volt), AP_data_cell$time, 1)/samplerate*1000
  AP_data_cell$cell<-rep(cellname[i],length(AP_data_cell$volt))
  AP_data_cell$genotype<-rep(curr_cell[2],length(AP_data_cell$volt))
  AP_data<-rbind(AP_data,AP_data_cell)
  #select indices that are above the firing threshold defined as 20 V/s
  AP_thres_ind<-which(AP_data_cell$first_deriv[1:601]>20)
  #calculate the difference between the indices >20 V/s and add one to the start
  AP_thres_diff<-c(2,diff(AP_thres_ind))
  #select the last index of the ones where thedifference to the previous index was >1
  AP_ind_thres<-AP_thres_ind[tail(which(AP_thres_diff>1),1)]
  fAHP_ind<-head(which(AP_data_cell$volt[600:(range_after+601)]==min(AP_data_cell$volt[600:(range_after+601)])),1)
  amplitude<-max(AP_data_cell$volt)-AP_data_cell$volt[fAHP_ind]
  max_ind<-head(which(AP_data_cell$volt==max(AP_data_cell$volt)),1)
  FWHA_ind<-which(AP_data_cell$volt>(amplitude/2+AP_data_cell$volt[fAHP_ind]))
  AP_properties_cell<-data.frame(
    "cell"=cellname[i],
    "genotype"=curr_cell[2],
    "rheobase"=AP_data_cell$curr[AP_ind_thres],
    "threshold"=AP_data_cell$volt[AP_ind_thres],
    "FWHA"=AP_data_cell$time[tail(FWHA_ind,1)]-AP_data_cell$time[head(FWHA_ind,1)],
    "fAHP"=AP_data_cell$volt[fAHP_ind]-mean(AP_data_cell$volt[1:100]),
    "risingetime"=AP_data_cell$time[max_ind]-AP_data_cell$time[AP_ind_thres],
    "repolarizingtime"=AP_data_cell$time[fAHP_ind]-AP_data_cell$time[max_ind],
    "amplitude"=amplitude
  )
  AP_properties<-rbind(AP_properties,AP_properties_cell)
  }
AP_data.summary <- AP_data %>%
  group_by(genotype,ind) %>%
  summarise(
    sd.volt = sd(volt, na.rm = TRUE),
    mn.volt = mean(volt),
    sd.firder=sd(first_deriv, na.rm = TRUE),
    mn.firder = mean(first_deriv),
    .groups="keep"
  )                  
  
  ggplot(AP_data,aes(volt,first_deriv, group=genotype, col=genotype))+
    geom_path()
    
p1<-  ggplot(data=AP_data,aes(volt,first_deriv))+
    geom_path(col="lightgrey")+
    geom_path(data=AP_data.summary,aes(mn.volt,mn.firder,col=genotype))+
    geom_errorbar(data=AP_data.summary,aes(mn.volt,mn.firder,xmin=mn.volt-sd.volt/sqrt(length(cellname)),xmax=mn.volt+sd.volt/sqrt(length(cellname)),ymin=mn.firder-sd.firder/sqrt(length(cellname)),ymax=mn.firder+sd.firder/sqrt(length(cellname)), col=genotype))+
    geom_errorbarh(data=AP_data.summary,aes(mn.volt,mn.firder,xmin=mn.volt-sd.volt/sqrt(length(cellname)),xmax=mn.volt+sd.volt/sqrt(length(cellname)),ymin=mn.firder-sd.firder/sqrt(length(cellname)),ymax=mn.firder+sd.firder/sqrt(length(cellname)), col=genotype))+
    scale_colour_manual(values = c("black", "blue","lightblue")) +
    ylab("first derivate [V/s]")+xlab("membran potential [mV]")+
    theme_classic()

ggsave(p1,width = 6, height = 6,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/phase_plot.png")
  
  ggplot(AP_data,aes(ind,volt, col=genotype))+
    geom_line(aes(group=cell))  
  ggplot(AP_data[AP_data$cell=="cell04",],aes(ind,volt, col=genotype))+
    geom_line(aes(group=cell)) 
  
p2<-ggplot(AP_properties,aes(genotype,rheobase,fill=genotype, col=genotype))+
    geom_boxplot()+
    geom_point()+
    scale_colour_manual(values = c("black", "blue","lightblue")) +
    scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    ylim(c(0,300))+
    ylab("rheobase [pA]")+
    theme_classic()
  wilcox.test(rheobase~genotype,AP_properties)
  ggsave(p2,width = 4, height = 4,
         file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/rheobase.png")
 
p3<-  ggplot(AP_properties,aes(genotype,threshold,fill=genotype, col=genotype))+
    geom_boxplot()+
    geom_point()+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    ylim(c(-60,0))+
    ylab("AP threshold [mV]")+
    theme_classic()
  wilcox.test(threshold~genotype,AP_properties)
  ggsave(p3,width = 4, height = 4,
         file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/threshold.png")
  
p4<-  ggplot(AP_properties,aes(genotype,FWHA,fill=genotype, col=genotype))+
    geom_boxplot()+
    geom_point()+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    ylab("FWHA [s]")+
    theme_classic()
  wilcox.test(FWHA~genotype,AP_properties)
  ggsave(p4,width = 4, height = 4,
         file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/FWHA.png")

p5<-  ggplot(AP_properties,aes(genotype,fAHP,fill=genotype, col=genotype))+
    geom_boxplot()+
    geom_point()+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    theme_classic()    
  wilcox.test(fAHP~genotype,AP_properties)
  ggsave(p5,width = 4, height = 4,
         file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/fAHP.png")
    
p6<-  ggplot(AP_properties,aes(genotype,risingetime,fill=genotype, col=genotype))+
    geom_boxplot()+
    geom_point()+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    theme_classic() 
  wilcox.test(risingetime~genotype,AP_properties)
  ggsave(p6,width = 4, height = 4,
         file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/risingetime.png")
    
p7<-  ggplot(AP_properties,aes(genotype,repolarizingtime,fill=genotype, col=genotype))+
    geom_boxplot()+
    geom_point()+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    theme_classic() 
wilcox.test(repolarizingtime~genotype,AP_properties,exact=T)  
ggsave(p7,width = 4, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/repolarizingtime.png")
  
p8<-  ggplot(AP_properties,aes(genotype,amplitude,fill=genotype, col=genotype))+
  geom_boxplot()+
  geom_point()+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  theme_classic() 
wilcox.test(amplitude~genotype,AP_properties,exact=T)  
ggsave(p8,width = 4, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/amplitude.png")
