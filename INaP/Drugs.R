library(tidyverse)
library(readxl)
library(xlsx)
library(cowplot)
library(viridis)
data<-read_xlsx("D:/Peter/Dokumente/Review_INaP/inap_ramp_pulse.xlsx")
data$Drug<-as.factor(data$Drug)
data.summary <- data %>%
  group_by(Drug) %>%
  summarise(
    sd.eff = sd(Effect, na.rm = TRUE),
    eff = mean(Effect),
    sd.dos=sd(Dose, na.rm = TRUE),
    dos = mean(Dose)
  )

data2<-read_xlsx("D:/Peter/Dokumente/Review_INaP/inap_slow.xlsx")
data2$Drug<-as.factor(data2$Drug)
data2.summary <- data2 %>%
  group_by(Drug) %>%
  summarise(
    sd.eff = sd(Effect, na.rm = TRUE),
    eff = mean(Effect),
    sd.dos=sd(Dose, na.rm = TRUE),
    dos = mean(Dose)
  ) 
data.summary$Drug<-factor(data.summary$Drug, levels = unique(c(data$Drug,data2$Drug)))
data2.summary$Drug<-factor(data2.summary$Drug, levels = unique(c(data$Drug,data2$Drug)))
palette<-viridis(22)
names(palette)=unique(c(data$Drug,data2$Drug))


a<-ggplot(data.summary, aes(dos, eff,xmin=dos-sd.dos,xmax=dos+sd.dos,ymin=eff-sd.eff,ymax=eff+sd.eff,color=Drug,shape=Drug))+
  scale_color_manual(values = palette,drop = FALSE)+
  scale_shape_manual(values = c(
    "Amiodarone" = 1,      
    "Cannabidiol"    =2, 
    "Carbamazepine"   =3,
    "Cenobamate"      =4,
    "Eslicarbazepine" =5,
    "Ethosuximide"    =6,
    "Gabapentin"      =7,
    "GS967"          =8,
    "Lacosamide"      =9,
    "Lamotrigene"     =10,
    "Lidocaine"       =11,
    "NBI-921352"      =12,
    "Oxcarbazepine"   =13,
    "Phenytoin"       =14,
    "PRAX-562"        =15,
    "Propofol"       =16,
    "Ranolazine"      =17,
    "Riluzole"        =18,
    "Rufinamide"      =19,
    "Topiramate"      =20,
    "Valproic acid"   =21,
    "Zonisamide"     =22),
    drop = FALSE
  )+
  geom_point(stroke=1.5)+
  guides(color = guide_legend(byrow = TRUE,ncol=1))+
  geom_errorbar()+
  geom_errorbarh()+
  theme_classic()+
  theme(text=element_text(size=10),
        legend.spacing.y= unit(-0.18,"cm"))+
  ylab(expression('Inhibition of I'[NaP]*'[%]'))+
  xlab(expression('log'[10]*' concentration [M]'))+
  ggtitle("a")+
  theme(plot.title = element_text(size=12, face ="bold"),plot.title.position = "plot")
#ggsave("D:/Peter/Dokumente/Review_INaP/ramp.svg",a)

      
b<-ggplot(data2.summary, aes(dos, eff,xmin=dos-sd.dos,xmax=dos+sd.dos,ymin=eff-sd.eff,ymax=eff+sd.eff,color=Drug,shape=Drug))+
  scale_color_manual(values = palette)+
  scale_shape_manual(values = c(
    "Amiodarone" = 1,      
    "Cannabidiol"    =2, 
    "Carbamazepine"   =3,
    "Cenobamate"      =4,
    "Eslicarbazepine" =5,
    "Ethosuximide"    =6,
    "Gabapentin"      =7,
    "GS967"          =8,
    "Lacosamide"      =9,
    "Lamotrigene"     =10,
    "Lidocaine"       =11,
    "NBI-921352"      =12,
    "Oxcarbazepine"   =13,
    "Phenytoin"       =14,
    "PRAX-562"        =15,
    "Propofol"       =16,
    "Ranozoline"      =17,
    "Riluzole"        =18,
    "Rufinamide"      =19,
    "Topiramate"      =20,
    "Valproic acid"   =21,
    "Ranezoline"      =22,
    "Zonisamide"     =23)
  )+
  geom_point(stroke=1.5)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  geom_errorbar()+
  geom_errorbarh()+
  theme_classic()+
  ylab(expression('shift of slow inactivation V'[0.5]*' [mV]'))+
  xlab(expression('log'[10]*' concentration [M]'))+
  scale_y_reverse()+
  theme(legend.position = "none",text=element_text(size=10))+
  ggtitle("b")+
  theme(plot.title = element_text(size=12, face ="bold"),plot.title.position = "plot")
#ggsave("D:/Peter/Dokumente/Review_INaP/slow.svg",b)

#ggarrange(a, b, ncol=2, common.legend = TRUE, legend="bottom")

legend <- cowplot::get_legend(a + theme(legend.position = "right"))
a<-a+ theme(legend.position = "none")
ab<-cowplot::plot_grid(a, legend,b, ncol = 3, rel_widths = c(0.8, 0.3, 0.6))
ggsave("D:/Peter/Dokumente/Review_INaP/figure3_20240411.pdf",ab, width = 17.4, height = 10, units = "cm")
