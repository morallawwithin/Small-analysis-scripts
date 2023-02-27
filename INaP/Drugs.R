library(tidyverse)
library(readxl)
data<-read_xlsx("D:/Peter/Dokumente/Review_INaP/inap.xlsx")
ggplt(data, aes(Dose, Effect, group=Drug))+
  