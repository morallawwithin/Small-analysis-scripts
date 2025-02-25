library(readABF)
library(tidyverse)
library(readxl)
library(ggprism)
library(pspline)
library(ggsignif)
library(ggbeeswarm)
library(data.table)
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys")

####################
##select the cells
####################
dataset<-"Cortex_L2&3_PN"#"CA1_PN"#"EC_L5PN"#"Cortex_L2&3_PN"
data <- read_excel(paste0(dataset,".xlsx"))
setwd(paste0("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/",dataset))
data<-data[data$protocol=="ramp",]
cells<-data[,c(2,3,5)]
cellname<-data$cell
##########
#prepare everything for the loop
##########
AP_data<-data.frame(matrix(ncol = 9, nrow = 0))
colnames(AP_data)<-c("volt","curr","time","ind","first_deriv","second_deriv","cell","genotype","age")
AP_properties<-data.frame(matrix(ncol = 10, nrow = 0))
colnames(AP_properties)<-c("cell","genotype","rheobase","threshold","FWHA","fAHP","risingetime","repolarizingtime","amplitude","age")
#Samplerate
samplerate<-100000 #/s

########
#loop starts here
#######
for ( i in 1:length(cellname)){
  #select the cell
  curr_cell<-cells[i,]
  #read data of one cell  
  data<-readABF(curr_cell$file)
  sweep_length<-length(data$data[[1]][,1])/samplerate
  command_curr  <-c(rep(0,0.07*samplerate),c(1:(samplerate*1.5))/(samplerate*1.5/200),rep(0,(0.5*samplerate-0.07*samplerate)))
  #plot(data$data[[1]][,2]*1000,type="l")
  #lines(command_curr,col="red")
  
  #find first AP
  AP_ind_raw<-head(which(data$data[[1]][,1]>0),1)
  range_after<-samplerate*0.03
  if (dataset=="CA1_PN"){
    range_after<-samplerate*0.013 
  }
  AP_data_cell<- data.frame(
    "volt" = c(data$data[[1]][(AP_ind_raw-samplerate*0.006):(AP_ind_raw+range_after),1]),
    "curr" = command_curr[(AP_ind_raw-samplerate*0.006):(AP_ind_raw+range_after)],
    "time" = (c(1:(sweep_length*samplerate))/samplerate)[(AP_ind_raw-samplerate*0.006):(AP_ind_raw+range_after)],
    "ind" = c(1:((samplerate*0.006+1)+range_after)))
  AP_data_cell$first_deriv<-predict(sm.spline(AP_data_cell$time, AP_data_cell$volt), AP_data_cell$time, 1)/1000
  AP_data_cell$second_deriv<-predict(sm.spline(AP_data_cell$time, AP_data_cell$first_deriv), AP_data_cell$time, 1)/1000
  AP_data_cell$cell<-rep(cellname[i],length(AP_data_cell$volt))
  AP_data_cell$genotype<-rep(curr_cell$genotype,length(AP_data_cell$volt))
  AP_data_cell$age<-rep(curr_cell$age,length(AP_data_cell$volt))
  AP_data<-rbind(AP_data,AP_data_cell)
  #select indices that are above the firing threshold defined as 20 V/s
  AP_thres_ind<-which(AP_data_cell$first_deriv[1:(samplerate*0.006+1)]>20)
  #calculate the difference between the indices >20 V/s and add one to the start
  AP_thres_diff<-c(2,diff(AP_thres_ind))
  #select the last index of the ones where thedifference to the previous index was >1
  AP_ind_thres<-AP_thres_ind[tail(which(AP_thres_diff>1),1)]
  fAHP_ind<-head(which(AP_data_cell$volt[(samplerate*0.006):(range_after+(samplerate*0.006+1))]==min(AP_data_cell$volt[(samplerate*0.006):(range_after+(samplerate*0.006+1))])),1)+(samplerate*0.006-1)
  amplitude<-max(AP_data_cell$volt)-AP_data_cell$volt[AP_ind_thres]
  max_ind<-head(which(AP_data_cell$volt==max(AP_data_cell$volt)),1)
  FWHA_ind<-which(AP_data_cell$volt>(amplitude/2+AP_data_cell$volt[fAHP_ind]))
  FWHA_diff<-diff(FWHA_ind)
  if(any(FWHA_diff>1)){
    FWHA_ind<-FWHA_ind[1:(which(FWHA_diff>1)-1)]
  }
  AP_properties_cell<-data.frame(
    "cell"=cellname[i],
    "genotype"=curr_cell$genotype,
    "rheobase"=AP_data_cell$curr[AP_ind_thres],
    "threshold"=AP_data_cell$volt[AP_ind_thres],
    "FWHA"=AP_data_cell$time[tail(FWHA_ind,1)]-AP_data_cell$time[head(FWHA_ind,1)],
    "fAHP"=AP_data_cell$volt[fAHP_ind]-mean(AP_data_cell$volt[1:(samplerate*0.001)]),
    "risingetime"=AP_data_cell$time[max_ind]-AP_data_cell$time[AP_ind_thres],
    "repolarizingtime"=AP_data_cell$time[fAHP_ind]-AP_data_cell$time[max_ind],
    "amplitude"=amplitude,
    "age"=curr_cell$age
  )
  AP_properties<-rbind(AP_properties,AP_properties_cell)
  }
setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/Cortex_L2&3_PN/P12-P16")
AP_data_all<-AP_data
AP_data<-AP_data_all[AP_data_all$age<17,]
AP_properties_all<-AP_properties
AP_properties<-AP_properties_all[AP_properties_all$age<17,]

#setwd("D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/Cortex_L2&3_PN/P17-P20")
#AP_data_all<-AP_data
#AP_data<-AP_data_all[AP_data_all$age>16,]
#AP_properties_all<-AP_properties
#AP_properties<-AP_properties_all[AP_properties_all$age>16,]

AP_data.summary <- AP_data %>%
  group_by(genotype,ind) %>%
  summarise(
    sd.volt = sd(volt, na.rm = TRUE),
    mn.volt = mean(volt),
    sd.firder=sd(first_deriv, na.rm = TRUE),
    mn.firder = mean(first_deriv),
    sd.secder=sd(second_deriv, na.rm = TRUE),
    mn.secder = mean(second_deriv),
    .groups="keep"
  )                  
  
#ggplot(AP_data,aes(volt,first_deriv, group=genotype, col=genotype))+
#    geom_path()
    
p1<-ggplot(data=AP_data,aes(volt,first_deriv))+
    geom_path(col="lightgrey")+
    geom_path(data=AP_data.summary,aes(mn.volt,mn.firder,col=genotype))+
    geom_errorbar(data=AP_data.summary,aes(mn.volt,mn.firder,xmin=mn.volt-sd.volt/sqrt(length(cellname)),xmax=mn.volt+sd.volt/sqrt(length(cellname)),ymin=mn.firder-sd.firder/sqrt(length(cellname)),ymax=mn.firder+sd.firder/sqrt(length(cellname)), col=genotype))+
    geom_errorbarh(data=AP_data.summary,aes(mn.volt,mn.firder,xmin=mn.volt-sd.volt/sqrt(length(cellname)),xmax=mn.volt+sd.volt/sqrt(length(cellname)),ymin=mn.firder-sd.firder/sqrt(length(cellname)),ymax=mn.firder+sd.firder/sqrt(length(cellname)), col=genotype))+
    scale_colour_manual(values = c("black", "blue","lightblue")) +
    ylab("first derivate [V/s]")+xlab("membran potential [mV]")+
    theme_prism(base_size = 14)
p1
ggsave(p1,width = 10, height = 6,
       file="phase_plot.png")
ggsave(p1,width = 10, height = 6,
       file="phase_plot.svg")


#second derivative  
p11<-ggplot(data=AP_data.summary,aes(ind/samplerate,mn.secder,col=genotype,group=genotype))+
  geom_path(data=AP_data,aes(ind/samplerate,second_deriv,group=cell),col="lightgrey")+
  geom_path(data=AP_data.summary,aes(ind/samplerate,mn.secder,col=genotype,group=genotype),size=1.5)+
  xlab("time [s]")+ylab("second derivative [V^2/s^2]")+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  xlim(0.0055,0.0065)+
  theme_prism(base_size = 14)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
p11
#mean volt, mean first derivative, mean second derivative
AP_data.summary2<-melt(setDT(AP_data.summary[,c(1,2,4,6,8)]), id.vars = c("genotype","ind"), variable.name = "mean")
p12<-ggplot(data=AP_data.summary2,aes(ind/samplerate,value))+
  geom_path(size=1.5,aes(color=genotype))+
  xlim(0.0055,0.008)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  facet_grid(rows = vars(mean), 
             #cols = vars(genotype),
             scales="free_y")+
  theme_prism()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")+
  xlab("time [s]")
p12
ggsave(p12,width = 4, height = 10,
       file="mean_AP_first_second.png")
ggsave(p12,width = 4, height = 10,
       file="mean_AP_first_second.svg")
#mean AP shape
p10<-ggplot(data=AP_data.summary,aes(ind/samplerate,mn.volt,col=genotype,group=genotype))+
  geom_path(data=AP_data,aes(ind/samplerate,volt,group=cell),col="lightgrey")+
  geom_path(data=AP_data.summary,aes(ind/samplerate,mn.volt,col=genotype,group=genotype),size=1.5)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  xlab("time [s]")+ylab("membran potential [mV]")+
  ylim(c(-80,50))+
  xlim(0,0.036)+
  geom_segment(aes(x = 0.022, y = -75, xend = 0.032, yend = -75),lwd=2,col= "black")+ #the length of the x bar is 1/8 of the total length of the recording
  geom_segment(aes(x = 0.032, y = -75, xend = 0.032, yend = -55),lwd=2,col= "black")   +           #the length of the y bar is 100 µV 
  theme_prism(base_size = 14)+
  geom_text(x=0.034, y=-65, label="20 mV",col= "black")+
  geom_text(x=0.028, y=-78, label="10 ms",col= "black")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

p10
ggsave(p10,width = 6, height = 10,
       file="mean_AP.png")
ggsave(p10,width = 6, height = 10,
       file="mean_AP.svg")
  # ggplot(AP_data,aes(ind,volt, col=genotype))+
  #   geom_line(aes(group=cell))  
  # ggplot(AP_data[AP_data$cell=="cell12",],aes(ind,volt))+ #05
  #   geom_line(aes(group=cell)) +
  #   geom_line(data=AP_data[AP_data$cell=="cell17",], aes(group=cell),col="blue") +#22
  #   theme_minimal()
  # #EC
  # ggplot(AP_data[AP_data$cell=="cell16",],aes(ind,volt))+ #05
  #   geom_line(aes(group=cell)) +
  #   geom_line(data=AP_data[AP_data$cell=="cell07",], aes(group=cell),col="blue") +#22
  #   theme_minimal()
  
p2<-ggplot(AP_properties,aes(genotype,rheobase,fill=genotype, col=genotype))+
    geom_boxplot()+
    geom_beeswarm(cex=5,size=4)+
    scale_colour_manual(values = c("black", "blue","lightblue")) +
    scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
    ylab("rheobase [pA]")+
    theme_prism(base_size = 14)+
    theme(legend.position = "none")
p2  
wilcox.test(rheobase~genotype,AP_properties)
  ggsave(p2,width = 3, height = 4,
         file="rheobase.png")
  ggsave(p2,width = 3, height = 4,
         file="rheobase.svg")
 
p3<-  ggplot(AP_properties,aes(genotype,threshold,fill=genotype, col=genotype))+
    geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    ylim(c(-60,0))+
    ylab("AP threshold [mV]")+
    theme_prism(base_size = 14)+
  theme(legend.position = "none")
p3
wilcox.test(threshold~genotype,AP_properties)
  ggsave(p3,width = 3, height = 4,
         file="threshold.png")
  ggsave(p3,width = 3, height = 4,
         file="threshold.svg")
  
p4<-  ggplot(AP_properties,aes(genotype,FWHA*1000,fill=genotype, col=genotype))+
    geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    ylab("FWHA [ms]")+
    theme_prism(base_size = 14)+
  theme(legend.position = "none")
p4  
wilcox.test(FWHA~genotype,AP_properties)
  ggsave(p4,width = 3, height = 4,
         file="FWHA.png")
  ggsave(p4,width = 3, height = 4,
         file="FWHA.svg")
p41<-ggplot(AP_properties_all,aes(age,FWHA*1000, col=as.factor(genotype),fill=as.factor(genotype)))+
    geom_point(shape=16, size=4)+
    geom_smooth(method='lm', formula= y~x)+
    scale_colour_manual(values = c("black", "blue")) +
    scale_fill_manual(values = c("lightgrey",rgb(191/255,191/255,1,1))) +
    theme_prism(base_size = 14)+
    xlim(c(12,20))+  
    xlab("age [d]") + ylab("FWHA [ms]")+
    theme(legend.position = "none") 
p41
fwha.mdl<-lm(FWHA~age*genotype,data=AP_properties)
anova(fwha.mdl)
ggsave(p41,width = 4, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/Cortex_L2&3_PN/fwha_vs_age.png")
p5<-  ggplot(AP_properties,aes(genotype,abs(fAHP),fill=genotype, col=genotype))+
    geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    theme_prism(base_size = 14)+
  ylab("afterhyperpolarization [mV]")+
  scale_y_continuous(expand = c(0, 0),limits = c(0,20))+
  theme(legend.position = "none")    
p5  
wilcox.test(fAHP~genotype,AP_properties)
  ggsave(p5,width = 3, height = 4,
         file="fAHP.png")
  ggsave(p5,width = 3, height = 4,
         file="fAHP.svg")
p55<-ggplot(AP_properties_all,aes(age,abs(fAHP), col=as.factor(genotype),fill=as.factor(genotype)))+
  geom_point(shape=16, size=4)+
  geom_smooth(method='lm', formula= y~x)+
  scale_colour_manual(values = c("black", "blue")) +
  scale_fill_manual(values = c("lightgrey",rgb(191/255,191/255,1,1))) +
  theme_prism(base_size = 14)+
  scale_y_continuous(expand = c(0, 0),limits = c(0,20))+
  xlim(c(12,20))+  
  xlab("age [d]") + ylab("afterhyperpolarization [mV]")+
  theme(legend.position = "none")  
p55      
ggsave(p55,width = 4, height = 4,
       file="D:/Peter/Analysis/KCNA2/P405L_Mice/E-Phys/Cortex_L2&3_PN/ahp_vs_age.png")
ahp.mdl<-lm(fAHP~age*genotype,data=AP_properties)
anova(ahp.mdl)
p6<-  ggplot(AP_properties,aes(genotype,risingetime,fill=genotype, col=genotype))+
    geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    theme_prism(base_size = 14)+
  theme(legend.position = "none") 
p6  
wilcox.test(risingetime~genotype,AP_properties)
  ggsave(p6,width = 3, height = 4,
         file="risingetime.png")
  ggsave(p6,width = 3, height = 4,
         file="risingetime.svg")  
p7<-  ggplot(AP_properties,aes(genotype,repolarizingtime,fill=genotype, col=genotype))+
    geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
    theme_prism(base_size = 14)+
  theme(legend.position = "none") 
p7
wilcox.test(repolarizingtime~genotype,AP_properties,exact=T)  
ggsave(p7,width = 3, height = 4,
       file="repolarizingtime.png")
ggsave(p7,width = 3, height = 4,
       file="repolarizingtime.svg")  
p8<-  ggplot(AP_properties,aes(genotype,amplitude,fill=genotype, col=genotype))+
  geom_boxplot()+
  geom_beeswarm(cex=5,size=4)+
  scale_colour_manual(values = c("black", "blue","lightblue")) +
  scale_fill_manual(values = c("white",rgb(191/255,191/255,1,1),"white"))+
  theme_prism(base_size = 14)+
  theme(legend.position = "none") 
p8
wilcox.test(amplitude~genotype,AP_properties,exact=T)  
ggsave(p8,width = 3, height = 4,
       file="amplitude.png")
ggsave(p8,width = 3, height = 4,
       file="amplitude.svg")
data1<-readABF(cells$file[1])
data1<-as.data.frame(data1)
data2<-readABF(cells$file[12])
data2<-as.data.frame(data2)
p9<-ggplot(data1,aes(`Time [s]`,`INcc 0 [mV]`))+
  geom_line()+
  geom_line(data=data2,color="blue")+
  ylim(c(-90,50))+
  theme_prism(base_size = 14)+
  xlab("time [s]") + ylab("memb. pot. [mV]")
ggsave(p9,
       file="examples_200pA_1s.png")  
