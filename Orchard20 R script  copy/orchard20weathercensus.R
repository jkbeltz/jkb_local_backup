weather=orchard20weather
census=finalcensuscountdata

library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(nlme)
library(lme4)

# extracting time
weather$time <- format(as.POSIXct(weather$TimeStamp ),format = "%H:%M:%S")

# extracting date
weather$date <- as.Date (weather$TimeStamp)

view(weather)
weather.sum=weather %>%
  group_by(date) %>%
  dplyr::summarise(N=n(),
            PTempC=mean(PTemp_C_Avg,na.rm = TRUE),
            PTempCsd=sd(PTemp_C_Avg,na.rm = TRUE),
            PtempCse= PTempCsd/sqrt(N),
            MaxAirTC=max(AirTC_Avg,na.rm = TRUE),
            MinAIRTC=min(AirTC_Avg,na.rm =TRUE),
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
xdata=unique(weather.sum)
View(xdata)


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
  
AirTCMin <-ggplot(wdata, aes(x=as.factor(date), y=MaxAirTC, group=1))+
  geom_line(color="Red")+
  #geom_linerange(aes(ymin=AirTC-AirTCse,ymax=AirTC+AirTCse))+
  ylab("Max/Min Daily Temperature C")+
  scale_x_discrete(name=" Experiment Date", 
                   breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
                   labels=c("Start","August", "September", "October", "November", "End")) +
  #xlab("Sample Timepoints")+
  #scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  #scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot()+
  theme(legend.position = "none")
AirTCMin

AirTCMinMax<-AirTCMin+geom_line(linetype= "solid" , data=wdata, aes(x=as.factor(date), y=MinAIRTC, group=1), color="Blue")
AirTCMinMaxCombined<-AirTCMinMax+geom_line(data=census.sum, aes(x=as.Date(Date), y=estavg, color=Cage_treatment))
AirTCMinMaxCombined  
  
  
  ######
  View(census)
census= census %>% 
  rename(
    pic.avg = `pic avg`)
    
  census$Date=as.Date(census$Date, format = "%m/%d/%y")
  
  census.sum=census %>%
    group_by(Date, Cage_treatment) %>%
    dplyr::summarise(N=n(),
              countavg=mean(Picture_total,na.rm = TRUE),
              countsd=sd(Picture_total,na.rm = TRUE),
              countse= countsd/sqrt(N),
              combestavg=mean(est_total,na.rm = TRUE),
              combestsd=sd(est_total,na.rm = TRUE),
              combestse=combestsd/sqrt(N^2), 
              picavgest=mean(est_avg,na.rm = TRUE),
              picavgsd=sd(est_avg,na.rm = TRUE),
              picavgse=picavgsd/sqrt(N),
              pictotavg=mean(est_tot,na.rm = TRUE),
              pictotsd=sd(est_tot,na.rm = TRUE),
              pictotse= pictotsd/sqrt(N),
              
              
    )
  ######
  View(censusnoATLB2)
  
  censusnoBLOOM=census[!(census$Cage_treatment == "Bloom"),]
  censusnoATLB<-census[(census$Cage_treatment == "Apple") | (census$Cage_treatment == "Bloom"),]
  
  censusnoATLB1=subset(censusnoATLB, censusnoATLB$Date == c("2020-07-21") | censusnoATLB$Date == c("2020-08-14"))
  censusnoATLB2=subset(censusnoATLB, censusnoATLB$Date == c("2020-09-25") | censusnoATLB$Date == c("2020-10-09") | censusnoATLB$Date == c("2020-10-18") | censusnoATLB$Date == c("2020-11-01") | censusnoATLB$Date == c("2020-11-10") | censusnoATLB$Date == c("2020-11-21"))
  

  censusAB<-lme(pic.avg ~ Date*Cage_treatment, random=~1|Cage_number/Date, data=censusnoATLB2)
  anova(censusAB)

  
  
  census$Date=as.Date(census$Date, format = "%m/%d/%y")
  
  cagecensus.sum=census %>%
    group_by(Date, Cage_number) %>%
    dplyr::summarise(N=n(),
              countavg=mean(Picture_total,na.rm = TRUE),
              countsd=sd(Picture_total,na.rm = TRUE),
              countse= countsd/sqrt(N),
              combestavg=mean(est_total,na.rm = TRUE),
              combestsd=sd(est_total,na.rm = TRUE),
              combestse=combestsd/sqrt(N), 
              picavgest=mean(est_avg,na.rm = TRUE),
              picavgsd=sd(est_avg,na.rm = TRUE),
              picavgse=picavgsd/sqrt(N),
              pictotavg=mean(est_tot,na.rm = TRUE),
              pictotsd=sd(est_tot,na.rm = TRUE),
              pictotse= pictotsd/sqrt(N),
              
              
    )
  #censusdata=unique(census.sum)
  View(census.sum)
  
  censusnoBLOOM.sum=census.sum[!(census.sum$Cage_treatment == "Bloom"),]
  censusnoATLB.sum<-census.sum[(census.sum$Cage_treatment == "Apple") | (census.sum$Cage_treatment == "Bloom"),]
  censusAB.sum<-census.sum[(census.sum$Cage_treatment == "AT") | (census.sum$Cage_treatment == "LB"),]
  
  
  censusbycage.sum=census %>%
    group_by(Date, Cage_treatment, Cage_number) %>%
    dplyr::summarise(N=n(),
              countavg=mean(Picture_total,na.rm = TRUE),
              estavg=mean(est_total,na.rm=TRUE), 
              picavgest=mean(est_avg,na.rm = TRUE),
              picavgsd=sd(est_avg,na.rm = TRUE),
              picavgse=picavgsd/sqrt(N)
              
    )
  
  censusdata=unique(census.sum)
  install.packages("lme4")
  library(lme4)
  library(nlme)
  ?lme
  
  
  censuscagenoBLOOM.sum=censusbycage.sum[!(censusbycage.sum$Cage_treatment == "Bloom"),]
  censuscagenoATLB.sum<-censusbycage.sum[(censusbycage.sum$Cage_treatment == "Apple") | (censusbycage.sum$Cage_treatment == "Bloom"),]

  censusAB<-lme(picavgest ~ Date*Cage_treatment, random=~1|Cage_number/Date, data=censuscagenoATLB.sum)
  anova(censusAB)
  
  censusATLBAP<-lme(picavgest ~ Date*Cage_treatment, random=~1|Cage_number/Date, data=censuscagenoBLOOM.sum)
  anova(censusATLBAP)
  
  
  
  library(cowplot)
  View(censusbycage.sum)
  #censusbycagenoATLB.sum=censusbycage.sum[!(censusbycage.sum$Cage_treatment = c("AT"))]
  #censusbycageAPBLOOM.sum<-censusbycage.sum[(censusbycage.sum$Cage_treatment == "Apple") | (censusbycage.sum$Cage_treatment == "Bloom"),]
  
  View(census.sum)
  
  CountAvgC <- ggplot(censuscagenoATLB.sum, aes(x=as.Date(Date), y=picavgest, by=Cage_number, color=Cage_treatment))+
    #geom_point(size=2, alpha=0.5)+
    geom_line(size=.5, alpha=0.65)+
    geom_linerange(aes(ymin=picavgest-picavgse,ymax=picavgest+picavgse), size=1)+
    ylab("Population Size")+
    xlab("Sample Date")+
    #scale_x_discrete(name=" Experiment Date", 
                     #breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
                     #labels=c("Start","August", "September", "October", "November", "End")) +
    scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT", "Bloom"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations", "Bloomington Populations"),
                      values = c("#848FA2","#058ED9","#CC2D35", "Black"))+
    xlab("Sample Timepoints")+
    #scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
    #scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
    theme_cowplot(12)+
    theme(legend.position = "none")
  
CountAvgC
  
combined= CountAvgC + geom_line(linetype= "solid" , data=dplrcombined, aes(x=as.Date(Date), y=MinAIRTC, group=1), color="Blue")
combined

#CountAvgC 
    #geom_line () 
    

a=geom_line(linetype= "solid" , data=wdata, aes(x=date, y=MinAIRTC, group=1), color="Blue")

a
  
  CountByCageAvgC <- ggplot(censusbycage.sum, aes(x=as.Date(Date), y=countavg, group= Cage_number, color=Cage_treatment))+
    geom_point(size=2)+
    geom_line(size=1)+
    #geom_linerange(aes(ymin=countavg-countse,ymax=countavg+countse))+
    ylab("Census Counts")+
    xlab("Sample Date")+
    scale_x_discrete(name=" Experiment Date", 
                     breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
                     labels=c("Start","August", "September", "October", "November", "End")) +
    scale_colour_manual(name = "Cage Treatment", 
                        breaks = c("Apple","LB", "AT", "Bloom"),
                        labels = c("Control Populations","LB+ Populations", "AT+ Populations", "Bloomington Populations"),
                        values = c("#848FA2","#058ED9","#CC2D35", "Black"))+
    #xlab("Sample Timepoints")+
    #scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
    #scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
    theme_cowplot()
    #theme(legend.position = "none")
  CountByCageAvgC

  
  
  
  
  library(dplyr)
  names(wdata)[1]= "Date"
  View(wdata)
  dplrcombinedtreatment <- left_join(wdata, census.sum, by="Date")
  dplrcombinedcage <- left_join(wdata, cagecensus.sum, by="Date")
  dplrcombinednoATLB <- left_join(wdata, censusnoATLB.sum, by="Date")
  dplrcombinednoBloom <- left_join(wdata, censusnoBLOOM.sum, by="Date")
  View(dplrcombined)
  
  AirTCMin <-ggplot(dplrcombined, aes(x=as.factor(Date), y=MaxAirTC, group=1))+
    geom_line(color="Red")+
    geom_linerange(aes(ymin=AirTC-AirTCse,ymax=AirTC+AirTCse))+
    ylab("Max/Min Daily Temperature C")+
    scale_x_discrete(name=" Experiment Date", 
                     breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
                     labels=c("Start","August", "September", "October", "November", "End")) +
    xlab("Sample Timepoints")+
    scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
    scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
    theme_cowplot()+
    theme(legend.position = "none")
  AirTCMin
  
  AirTCMinMax<-AirTCMin+geom_line(linetype= "solid" , data=dplrcombined, aes(x=as.factor(Date), y=MinAIRTC, group=1), color="Blue")
  AirTCMinMax ###### NEED ANOTHER AXIS 
  
  
  
  
  coeff <- 5000
  coeff2 <-100
  
  # A few constants
  temperatureColor <- "#69b3a2"
  priceColor <- rgb(0.2, 0.6, 0.9, 1)
  
  test=ggplot(dplrcombined, aes(x=as.factor(Date))) +
    geom_line( aes(y=MinAIRTC), size=1, group=1, color="Blue", alpha=.5) + 
    geom_line( aes(y=MaxAirTC), size=1, group=1, color="Red", alpha=.5) +
    geom_line( aes(y=combestavg / coeff), size=2, group=dplrcombined$Cage_treatment,  ) +
    scale_x_discrete(
      name=" Experiment Date", 
      breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
      labels=c("Start","August", "September", "October", "November", "End"))+
    
    scale_y_continuous(
      # Features of the first axis
      name = "Min/ Max Daily Temperature (Celsius °)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coeff, name="Population Size"))+
    theme_cowplot()
test
    
 #   theme(
  #    axis.title.y = element_text(color = temperatureColor, size=13),
  #    axis.title.y.right = element_text(color = priceColor, size=13)
  #  ) +
  #  
  #  ggtitle("Temperature and Population")







values = c("#848FA2","#CC2D35", "Black","#058ED9")
group.colors <- c(Apple = "#848FA2", AT = "#CC2D35", Bloom ="Black", LB = "#058ED9")

view(dplrcombined)
library(dplyr)
dplrcombined = dplrcombined %>%
  mutate(color = case_when(
    endsWith(Cage_treatment, "e") ~ "#848FA2",
    endsWith(Cage_treatment, "T") ~ "#CC2D35",
    endsWith(Cage_treatment, "m") ~ "Black",
    endsWith(Cage_treatment, "B") ~ "#058ED9"
  ))

dplrcombinedAPBLOOM<-dplrcombined[(dplrcombined$Cage_treatment == "Apple") | (dplrcombined$Cage_treatment == "Bloom"),]

test2=ggplot(dplrcombinednoBloom, aes(x=as.factor(Date))) +
  geom_line(aes(y=combestavg, group=dplrcombinednoBloom$Cage_treatment, color=dplrcombinednoBloom$Cage_treatment), size=1.3)+
  geom_errorbar(aes(ymin=combestavg-combestse,ymax=combestavg+combestse, color=dplrcombinednoBloom$Cage_treatment), size=.5)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT", "Bloom"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations", "Bloomington Populations"),
                      values = c("#848FA2","#058ED9","#CC2D35", "Black"))+
  #geom_linerange(aes(ymin=combestavg-combestse,ymax=combestavg+combestse), size=.5)+
  geom_line( aes(y=MinAIRTC * coeff), size=2, group=1, color="blue", alpha=.3) + 
  geom_line( aes(y=MaxAirTC * coeff), size=2, group=1, color="orange", alpha=.3) +
  scale_x_discrete(
    name=" Experiment Date", 
    breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
    labels=c("Start","August", "September", "October", "November", "End"))+
  
 scale_y_continuous(
    # Features of the first axis
    name = "Population Size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Min/ Max Daily Temperature (C°)"))+
  theme_cowplot()

test2+ theme(legend.position='none')
##### 
view(dplrcombinednoATLB)



test3=ggplot(dplrcombinednoATLB, aes(x=as.factor(Date))) +
  geom_line(aes(y=combestavg, group=dplrcombinednoATLB$Cage_treatment, color=dplrcombinednoATLB$Cage_treatment), size=1.3)+
  geom_errorbar(aes(ymin=combestavg-combestse,ymax=combestavg+combestse, color=dplrcombinednoATLB$Cage_treatment), size=.5)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT", "Bloom"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations", "Bloomington Populations"),
                      values = c("red","#058ED9","#CC2D35", "black"))+
  #geom_linerange(aes(ymin=combestavg-combestse,ymax=combestavg+combestse), size=.5)+
  geom_line( aes(y=MinAIRTC * coeff), size=2, group=1, color="blue", alpha=.3) + 
  geom_line( aes(y=MaxAirTC * coeff), size=2, group=1, color="orange", alpha=.3) +
  scale_x_discrete(
    name=" Experiment Date", 
    breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
    labels=c("Start","August", "September", "October", "November", "End"))+
  
  scale_y_continuous(
    # Features of the first axis
    name = "Population Size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Min/ Max Daily Temperature (C°)"))+
  theme_cowplot()

test3+ theme(legend.position='none')



######by count avg
dplrcombinedAPBLOOM<-dplrcombined[(dplrcombined$Cage_treatment == "Apple") | (dplrcombined$Cage_treatment == "Bloom"),]

test2=ggplot(dplrcombinednoBloom, aes(x=as.factor(Date))) +
  geom_line(aes(y=combestavg, group=dplrcombinednoBloom$Cage_treatment, color=dplrcombinednoBloom$Cage_treatment), size=1.3)+
  geom_errorbar(aes(ymin=combestavg-combestse,ymax=combestavg+combestse, color=dplrcombinednoBloom$Cage_treatment), size=.5)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT", "Bloom"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations", "Bloomington Populations"),
                      values = c("#848FA2","#058ED9","#CC2D35", "Black"))+
  #geom_linerange(aes(ymin=combestavg-combestse,ymax=combestavg+combestse), size=.5)+
  geom_line( aes(y=MinAIRTC * coeff2), size=2, group=1, color="blue", alpha=.3) + 
  geom_line( aes(y=MaxAirTC * coeff2), size=2, group=1, color="orange", alpha=.3) +
  scale_x_discrete(
    name=" Experiment Date", 
    breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
    labels=c("Start","August", "September", "October", "November", "End"))+
  
  scale_y_continuous(
    # Features of the first axis
    name = "Population Size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Min/ Max Daily Temperature (C°)"))+
  theme_cowplot()

test2+ theme(legend.position='none')
##### 
view(dplrcombinednoATLB)



test3=ggplot(dplrcombinednoATLB, aes(x=as.factor(Date))) +
  geom_line(aes(y=picavgest, group=dplrcombinednoATLB$Cage_treatment, color=dplrcombinednoATLB$Cage_treatment), size=1.3)+
  geom_errorbar(aes(ymin=picavgest-picavgse,ymax=picavgest+picavgse, color=dplrcombinednoATLB$Cage_treatment), size=.5)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT", "Bloom"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations", "Bloomington Populations"),
                      values = c("red","#058ED9","#CC2D35", "black"))+
  #geom_linerange(aes(ymin=combestavg-combestse,ymax=combestavg+combestse), size=.5)+
  geom_line( aes(y=MinAIRTC * coeff), size=2, group=1, color="blue", alpha=.3) + 
  geom_line( aes(y=MaxAirTC * coeff), size=2, group=1, color="orange", alpha=.3) +
  scale_x_discrete(
    name=" Experiment Date", 
    breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
    labels=c("Start","August", "September", "October", "November", "End"))+
  
  scale_y_continuous(
    # Features of the first axis
    name = "Population Size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Min/ Max Daily Temperature (C°)"))+
  theme_cowplot()

test3+ theme(legend.position='none')

