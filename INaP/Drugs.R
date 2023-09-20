library(tidyverse)
library(readxl)
library(xlsx)
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


a<-ggplot(data.summary, aes(dos, eff,xmin=dos-sd.dos,xmax=dos+sd.dos,ymin=eff-sd.eff,ymax=eff+sd.eff,color=Drug,shape=Drug))+
  scale_shape_manual(values=1:nlevels(data$Drug))+
  geom_point(aes(size=20))+
  geom_errorbar()+
  geom_errorbarh()+
  theme_classic()+
  ylab(expression('Inhibition of I'[NaP]*'[%]'))+
  xlab(expression('log'[10]*' concentration [M]'))+
  theme(legend.position = "none")
ggsave("D:/Peter/Dokumente/Review_INaP/ramp.svg",a)

data2<-read_xlsx("D:/Peter/Dokumente/Review_INaP/inap.xlsx")
data2$Drug<-as.factor(data2$Drug)
data2.summary <- data2 %>%
  group_by(Drug) %>%
  summarise(
    sd.eff = sd(Effect, na.rm = TRUE),
    eff = mean(Effect),
    sd.dos=sd(Dose, na.rm = TRUE),
    dos = mean(Dose)
  )       
b<-ggplot(data2.summary, aes(dos, eff,xmin=dos-sd.dos,xmax=dos+sd.dos,ymin=eff-sd.eff,ymax=eff+sd.eff,color=Drug,shape=Drug))+
  scale_shape_manual(values=1:nlevels(data$Drug))+
  geom_point(aes(size=20))+
  geom_errorbar()+
  geom_errorbarh()+
  theme_classic()+
  ylab(expression('shift of slow inactivation V'[0.5]*' [mV]'))+
  xlab(expression('log'[10]*' concentration [M]'))+
  scale_y_reverse()
ggsave("D:/Peter/Dokumente/Review_INaP/slow.svg",b)
