weather=orchard20weather

library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)

# extracting time
weather$time <- format(as.POSIXct(
  weather$TimeStamp ),format = "%H:%M:%S")

# extracting date
weather$date <- as.Date (weather$TimeStamp)

weather.sum=weather %>%
  group_by(date) %>%
  summarise(N=n(),
            PTempC=mean(PTemp_C_Avg,na.rm = TRUE),
            PTempCsd=sd(PTemp_C_Avg,na.rm = TRUE),
            PtempCse= PTempCsd/sqrt(N),
            AirTC=mean(AirTC_Avg,na.rm = TRUE),
            AirTCsd=sd(AirTC_Avg,na.rm = TRUE),
            AirTCse= AirTCsd/sqrt(N),
            RH=mean(`RH (percent)`,na.rm = TRUE),
            RHsd=sd(`RH (percent)`,na.rm = TRUE),
            RHse=RHsd/sqrt(N),
            Temp_C=mean(Temp_C,na.rm = TRUE),
            Temp_Csd=sd(Temp_C,na.rm = TRUE),
            Temp_Cse= Temp_Csd/sqrt(N),
            Temp_Cavg=mean(Temp_C_Avg,na.rm = TRUE),
            Temp_Cavgsd=sd(Temp_C_Avg,na.rm = TRUE),
            Temp_Cavgse= Temp_Cavgsd/sqrt(N),
            PARden=mean(PAR_Den_Avg,na.rm = TRUE),
            PARdensd=sd(PAR_Den_Avg,na.rm = TRUE),
            PARdense= PARdensd/sqrt(N),
            PARtot=mean(PAR_Tot_Tot,na.rm = TRUE),
            PARtotsd=sd(PAR_Tot_Tot,na.rm = TRUE),
            PARtotse= PARtotsd/sqrt(N),
            WS=mean(WS_ms_Avg,na.rm =TRUE),
            WSsd=sd(WS_ms_Avg,na.rm = TRUE),
            WSse= WSsd/sqrt(N),
            Rain=mean(Rain_mm_Tot,na.rm = TRUE),
            Rainsd=sd(Rain_mm_Tot,na.rm = TRUE),
            Rainse=Rainsd/sqrt(N)
  )
wdata=unique(weather.sum)
View(wdata)


AirTCw <- ggplot(wdata, aes(x=as.factor(date), y=AirTC, group=1))+
  geom_point()+
  geom_line()+
  geom_linerange(aes(ymin=AirTC-AirTCse,ymax=AirTC+AirTCse))+
  ylab("Average Daily Temperature")+
  #xlab("Sample Timepoints")+
  #scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  #scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot()+
  theme(legend.position = "none")
  AirTCw
  