library (plyr)
library(tidyverse)
library(readr)
library(visreg)
library(egg)
library(lme4)
library(multcomp)


setwd("D:/Peter/Analysis/CACNA1A") #adjust working directory
data <- read.csv("Phenomaster_G353R_20221203_corr.csv", header=TRUE)#adjust file
data$Box_date <- paste0(as.character(data$Box), "_221203")#adjust date


#data <- rbind(data1,data2)

data$Box<-as.factor(data$Box)

#delete empty boxes
#data<-data[!(data$Box_date %in% c("3_220923","3_220927","4__220513","5_220513","1_220513","1_220401","2_220401","6_220401")),]

#set groups
#Group1 <- c("1_220923","2_220923","4_220923","5_220923","6_220923","7_220923")#adjust your groups
#group_bin<-as.integer(data$Box_date %in% Group1)+1
#data$group<-as.factor(group_bin)
#levels(data$group)[levels(data$group)=="1"] <- "Without 4-AP"
#levels(data$group)[levels(data$group)=="2"] <- "4-AP"

#set genotype
wildytype <- c("2_221203","5_221203","7_221203")#adjust your groups
genotype_bin<-as.integer(data$Box_date %in% wildytype)+1
data$genotype<-as.factor(genotype_bin)
levels(data$genotype)[levels(data$genotype)=="1"] <- "P353R"
levels(data$genotype)[levels(data$genotype)=="2"] <- "wt"

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
gene = as.factor(genotype_bin)
)

levels(bin$gene)[levels(bin$gene)=="1"] <- "P353R"
levels(bin$gene)[levels(bin$gene)=="2"] <- "wt"
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
#act.mdl<-lm(median_activity ~ group*cos(2*pi*time_mdl/(24*60/interval))+I(time_mdl^2),data=bin)#this is a linear model with a poly fit for time

act.mdl<-lm(median_activity ~ gene*cos(2*pi*time_mdl/(24*60/interval)),data=bin)#this is a linear model with a cos fit for time
#summary(glht(act.mdl, linfct=mcp(gene="Tukey")))
summary(act.mdl)
vis<-visreg(act.mdl,"time_mdl",by="gene",overlay=TRUE,ylim=c(0,800))
a<-ggplot(vis$fit, aes(time_mdl, visregFit))+
  geom_boxplot(data=filter(vis$res, gene == "P353R"), aes(time_mdl, visregRes,group=as.factor(time_mdl),fill="P353R"),position = position_nudge(x=-0.1),outlier.shape = NA)+
  geom_boxplot(data=filter(vis$res, gene == "wt"), aes(time_mdl, visregRes,group=as.factor(time_mdl),fill="wt"),position = position_nudge(x=0.1),outlier.shape = NA)+
  scale_fill_manual("", 
                    breaks = c("P353R", "wt"),
                    values = c("deepskyblue", "grey")) +
  geom_line(data=filter(vis$fit, gene == "P353R"),colour='blue', size=1)+
  geom_ribbon(data=filter(vis$fit, gene == "P353R"),aes(ymin=visregLwr, ymax=visregUpr), fill='blue',alpha=.3)+
  geom_line(data=filter(vis$fit, gene == "wt"),colour='black', size=1)+
  geom_ribbon(data=filter(vis$fit, gene == "wt"),aes(ymin=visregLwr, ymax=visregUpr), fill='darkgrey',alpha=.3)+
  labs(y = "Median Velocity by h [cm/min]")+
  labs(x="Time of Day")+
  scale_x_continuous(breaks=c(4,8,12,16,20,24), labels=c("4:00","8:00","12:00","16:00","20:00","24:00"),limits=c(0.5,24.5))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  theme(legend.position=c(.9,.85))

ggsave(a,file="Median_velocity_by_h.png")

act.mdl2<-lm(sum_activity ~ gene*cos(2*pi*time_mdl/(24*60/interval)),data=bin)#this is a linear model with a cos fit for time
#summary(glht(act.mdl, linfct=mcp(gene="Tukey")))
summary(act.mdl2)
vis2<-visreg(act.mdl2,"time_mdl",by="gene",overlay=TRUE)
aa<-ggplot(vis2$fit, aes(time_mdl, visregFit))+
  geom_boxplot(data=filter(vis2$res, gene == "P353R"), aes(time_mdl, visregRes,group=as.factor(time_mdl),fill="P353R"),position = position_nudge(x=-0.1),outlier.shape = NA)+
  geom_boxplot(data=filter(vis2$res, gene == "wt"), aes(time_mdl, visregRes,group=as.factor(time_mdl),fill="wt"),position = position_nudge(x=0.1),outlier.shape = NA)+
  scale_fill_manual("", 
                      breaks = c("P353R", "wt"),
                      values = c("deepskyblue", "grey")) +
  geom_line(data=filter(vis2$fit, gene == "P353R"),colour='blue', size=1)+
  geom_ribbon(data=filter(vis2$fit, gene == "P353R"),aes(ymin=visregLwr, ymax=visregUpr), fill='blue',alpha=.3)+
  geom_line(data=filter(vis2$fit, gene == "wt"),colour='black', size=1)+
  geom_ribbon(data=filter(vis2$fit, gene == "wt"),aes(ymin=visregLwr, ymax=visregUpr), fill='darkgrey',alpha=.3)+
  labs(y = "Total distance per h [cm]")+
  labs(x="Time of Day")+
  scale_x_continuous(breaks=c(4,8,12,16,20,24), labels=c("4:00","8:00","12:00","16:00","20:00","24:00"),limits=c(0.5,24.5))+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position=c(.9,.75))
ggsave(aa,file="Total_distance_per_h.png")

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
  box=unique(data$Box),
  gene= as.factor(as.integer(data$Box_date %in% wildytype)+1))
levels(data.metabolics$gene)[levels(data.metabolics$gene)=="1"] <- "P353R"
levels(data.metabolics$gene)[levels(data.metabolics$gene)=="2"] <- "wt"

feed<-t.test(feed~gene,data=data.metabolics)
group.col=c("deepskyblue","sienna1")
d<-ggplot(data.metabolics,aes(x=gene,y=feed, fill=gene))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=group.col)+
  geom_point(position=position_jitter(0.1))+
  labs(y = "Food Intake per day [g]")+
  labs(x="")+
  theme(legend.position = "none")
drink<-t.test(drink~gene,data=data.metabolics)
f<-ggplot(data.metabolics,aes(x=gene,y=drink, fill=gene))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=group.col)+
  geom_point(position=position_jitter(0.1))+
  labs(y = "Water Intake per day [ml]")+
  labs(x="")+
  theme(legend.position = "none")
b<-ggarrange(f,d,ncol=1)
ggarrange(a,b,ncol=2,widths=c(0.8,0.2))
