library(tidyverse)
library(readr)
library (plyr)
library(visreg)


setwd("D:/Peter/Analysis/KCNA2/P405L_Mice") #adjust working directory
data1 <- read.csv("Phenomaster_20220401.csv", header=TRUE)#adjust file
data1$Box <- paste0(as.character(data1$Box), "_220401")#adjust date
data2 <- read.csv("Phenomaster_20220404.csv", header=TRUE)
data2$Box <- paste0(as.character(data2$Box), "_220404")
data3 <- read.csv("Phenomaster_20220408.csv", header=TRUE)
data3$Box <- paste0(as.character(data3$Box), "_220408")
data4 <- read.csv("Phenomaster_20220513.csv", header=TRUE)
data4$Box <- paste0(as.character(data4$Box), "_220513")
data5 <- read.csv("Phenomaster_20220516.csv", header=TRUE)
data5$Box <- paste0(as.character(data5$Box), "_220516")
data <- rbind(data1,data2)
data <- rbind(data,data3)#adjust if you have more experiments
data <- rbind(data,data4)#adjust if you have more experiments
data <- rbind(data,data5)#adjust if you have more experiments
data$Box<-as.factor(data$Box)

#delete empty boxes
data<-data[!(data$Box %in% c("5_220408","2_220516")),]

#set groups
Group1 <- c("3_220401","4_220401","5_220401","7_220401","2_220404","3_220404","6_220404","1_220408","4_220408","2_220513","3_220513","6_220513","7_220513","3_220516","5_220516","6_220516")#adjust your groups
group_bin<-as.integer(data$Box %in% Group1)+1
data$group<-as.factor(group_bin)

#set Up the data container (bin)
n_min = 1
diff_activity_bin <- data$DistD
time_bin <- data$Time
data$time<-as.POSIXlt(time_bin, tz="", tryFormats = "%H:%M:%OS" )
name_bin <-data$Box
bin <- data.frame(
diff_activity = diff_activity_bin,
time =as.POSIXlt(time_bin, tz="", tryFormats = "%H:%M:%OS" ),
time_disp=time_bin,
name =name_bin,
group = as.factor(group_bin)
)

#generate timestamps every interval min
start <- as.POSIXct("00:00:00",tryFormats = "%H:%M:%OS")
interval <- 60
end <- start + as.difftime(1, units="days")
timepoints <- seq(from=start, by=interval*60, to=end)
time_logical <- bin$time %in% timepoints       
#sum acitivy between timestamps
sum_acitivity_bin<- rep(NaN,length(time_bin))#this creates a matrix with empty values
sum_acitivity <- tapply(data$DistD, cumsum(time_logical), sum)#this aplies a function (in this case sum) to all timepoints between chosen intervals
sum_acitivity_bin[time_logical]<-sum_acitivity[-1]
bin$sum_activity <- sum_acitivity_bin
#disp median activity between timestamps
median_acitivity_bin<- rep(NaN,length(time_bin))
median_acitivity_var <- tapply(data$DistD, cumsum(time_logical), median)
median_acitivity_bin[time_logical]<-median_acitivity_var[-1]
bin$median_activity <- median_acitivity_bin
#any_activity between timestamps
any_acitivity_bin<- rep(NaN,length(time_bin))
any_acitivity <- tapply(data$DistD, cumsum(time_logical), function(c)sum(c!=0)/interval)
any_acitivity_bin[time_logical]<-any_acitivity[-1]
bin$any_activity <- any_acitivity_bin

#plot data at timestamps
bin <-bin[time_logical,]#this deletes all non interval timepoints (nothing of value is lost)
bin$time_disp <- as.factor (bin$time_disp)#need this for plotting and testing
bin$time_mdl<-cos((as.integer(bin$time_disp)-8)/12*pi)
bin$time_mdl<-as.integer(bin$time_disp)
#Exclusion of the first 8 intervals,just delete this section if you want to analyze the full data
exclusion<-split(seq_along(bin$name),bin$name)
exclusion_indices<-sapply(exclusion, "[", 1:8) #change 8 to whatever you need
exclusion_logical<-!(seq(bin$name) %in% exclusion_indices) #if you only want the first 8 intervals remove the !
bin<-bin[exclusion_logical,]


#median_activity


##Statistics
#lets first look at the distribution of qantiles
act.mdl<-lm(median_activity ~ group*(I(time_mdl^4)+I(time_mdl^2)),data=bin)#this is a linear model with a poly fit for time
act.mdl<-lm(median_activity ~ group*cos(2*pi*time_mdl/24),data=bin)#this is a linear model with a poly fit for time
summary(act.mdl)
vis<-visreg(act.mdl,"time_mdl",by="group",overlay=TRUE,ylim=c(0,400))
##Plot
a<-ggplot(filter(vis$fit, group == 1), aes(time_mdl, visregFit))+
  geom_boxplot(data=filter(vis$res, group == 1), aes(time_mdl, visregRes,group=as.factor(time_mdl)),position = position_nudge(x=-0.1),fill="deepskyblue",outlier.shape = NA)+
  geom_boxplot(data=filter(vis$res, group == 2), aes(time_mdl, visregRes,group=as.factor(time_mdl)),position = position_nudge(x=0.1),fill="sienna1",outlier.shape = NA)+
  #geom_point(data=filter(vis$res, group == 1), aes(time_mdl, visregRes,group=as.factor(time_mdl)), size=0.5, alpha=.3, position=position_jitter(0.4), colour='blue')+
  #geom_point(data=filter(vis$res, group == 2), aes(time_mdl, visregRes,group=as.factor(time_mdl)), size=0.5, alpha=.3, position=position_jitter(0.4), colour='red')+
  geom_line(colour='blue', size=1)+
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill='blue',alpha=.3)+
  geom_line(data=filter(vis$fit, group == 2),colour='red', size=1)+
  geom_ribbon(data=filter(vis$fit, group == 2),aes(ymin=visregLwr, ymax=visregUpr), fill='red',alpha=.3)+
  ylim(-10,200)

#Food and Drink

data.drink<-lapply(as.list(unique(data$Box)), function(x){
  max(data$Drink[data$Box==x])}  )
data.drink<-unlist(data.drink)
data.feed<-lapply(as.list(unique(data$Box)), function(x){
  max(data$Feed[data$Box==x])}  )  
data.feed<-unlist(data.feed)
data.metabolics<-data.frame(
  drink=data.drink/71*24,
  feed=data.feed/71*24,
  group=as.factor(as.integer(unique(data$Box) %in% Group1)+1),
  box=unique(data$Box))
data.metabolics$drink[[c(11)]]<-NA #remove false drink values due to bottle emptying on the floor 
data.metabolics$drink[[c(16)]]<-NA 
feed<-t.test(feed~group,data=data.metabolics)
group.col=c("deepskyblue","sienna1")
d<-ggplot(data.metabolics,aes(x=group,y=feed, fill=group))+
  geom_boxplot()+
  scale_fill_manual(values=group.col)+
  geom_point(position=position_jitter(0.1))+
  labs(y = "Food Intake per day [g]")+
  labs(x="")+
  theme(legend.position = "none")
drink<-t.test(drink~group,data=data.metabolics)
f<-ggplot(data.metabolics,aes(x=group,y=drink, fill=group))+
  geom_boxplot()+
  scale_fill_manual(values=group.col)+
  geom_point(position=position_jitter(0.1))+
  labs(y = "Water Intake per day [ml]")+
  labs(x="")+
  theme(legend.position = "none")
