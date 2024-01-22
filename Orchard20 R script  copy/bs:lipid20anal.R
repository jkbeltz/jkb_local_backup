library(vegan)
library(ecodist)
library(reshape2)
library(ggplot2)
library(gdata)
library(plyr)
library(compute.es)
library(lme4)
library(nlme)
library(dplyr)
library(car)
#install.packages("MASS", "cowplot")
#install.packages("cowplot")
library(MASS)
library(ggplot2)
library(cowplot)
library(tidyr)
library(multcomp)

bslipid20=bs_lipid20anal_Sheet7
bslipid20$arsqperclipid=asin(sqrt(bslipid20$'%lipid'/100))
View(bslipid20)

Founder=bslipid20[(bslipid20$Pop=="Founder"),]
names(Founder)[2]='timepoint'
names(Founder)[5]='cagetreatment'
names(Founder)[7]='dryweight'
names(Founder)[6]='phenotreatment'
names(Founder)[8]='lipidweight'
names(Founder)[9]='perclipid'
View(Founder)

Founder.DWLP.sum=Founder %>%
  group_by(cagetreatment) %>%
  #filter(!(Timepoint > 1 )) %>% #Removing apple treatments from bloomington cages
  #filter(!(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing bloomington treatments from apple cages
  #filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
 #filter(!(Timepoint== 1))%>% ## remove timepoint 1
  summarise(N=n(),
            treatment=as.factor(cagetreatment),
            DW=mean(dryweight,na.rm = TRUE),
            DWsd=sd(dryweight,na.rm = TRUE),
            DWse= DWsd/sqrt(N),
            ASPL=mean(arsqperclipid,na.rm = TRUE),
            ASPLsd=sd(arsqperclipid,na.rm = TRUE),
            ASPLse= ASPLsd/sqrt(N),
  )
founddata=unique(Founder.DWLP.sum)
View(founddata)
#### FOUNDER GEOMPOINT DFs Aggregate ####
foundermean <-aggregate(Founder[c("dryweight","lipidweight", "arsqperclipid")], by=Founder[c("cagetreatment","timepoint", "phenotreatment" )], FUN=mean, na.rm=FALSE)
foundermean2 <-aggregate(Founder[c("dryweight","lipidweight", "arsqperclipid")], by=Founder[c("cagetreatment","timepoint" )], FUN=mean, na.rm=TRUE)
Applebloomfoundermean=foundermean2[foundermean2$cagetreatment %in% c("Apple", "Bloom"), ]  
#View(Applebloomfoundermean)
AppleATLBfoundermean=foundermean2[foundermean2$cagetreatment %in% c("Apple", "AT", "LB"), ]
AppleBloomATLBfoundermean=foundermean2[foundermean2$cagetreatment %in% c("Apple","Bloom", "AT", "LB"), ]
View(AppleATLBfoundermean)
Applefoundermean=foundermean2[foundermean2$cagetreatment %in% c("Apple"), ]   
Bloomfoundermean=foundermean2[foundermean2$cagetreatment %in% c( "Bloom"), ]  
ATaxfoundermean=foundermean[foundermean$cagetreatment %in% c("Axenic", "AT"), ]   
#View(ATaxfoundermean)
LBaxfoundermean=foundermean[foundermean$cagetreatment %in% c( "Axenic", "LB"), ]  



Evo=bslipid20[(bslipid20$Pop=="Evo"),]

Evo=as.data.frame(Evo)
names(Evo)[2]='timepoint'
names(Evo)[5]='cagetreatment'
names(Evo)[6]='phenotreatment'
names(Evo)[7]='dryweight'
names(Evo)[8]='lipidweight'
names(Evo)[9]='perclipid'

Evo$`Cage Treatment`<-factor(Evo$cagetreatment)
Evo$`Pheno-treatment`<-factor(Evo$phenotreatment)
Evo$`Time Point`<-factor(Evo$timepoint)

Evo=Evo[!(Evo$timepoint == "TP5"),]
Evol=Evo[!(Evo$cagetreatment == "Bloom"),]
Evol=Evol[!(Evol$phenotreatment == "Bloom"),]
Evol=Evol[!(Evol$phenotreatment == "Axenic"),]
View(Evol)

install.packages("lme4")

install.packages("nlme")
library(nlme)
library(lme4)
DW1<-lme(dryweight ~ timepoint*cagetreatment, random=~1|Cage/timepoint, data=Evol)
anova(DW1) 

LW1<-lme(lipidweight ~ timepoint*cagetreatment, random=~1|Cage/timepoint, data=Evol)
anova(LW1) 

PL1<-lme(perclipid ~ timepoint*cagetreatment, random=~1|Cage/timepoint, data=Evol)
anova(PL1) 


bloomevo=Evo[(Evo$cagetreatment =="Bloom"),]
bloomevomean <-aggregate(bloomevo[c("dryweight","lipidweight", "arsqperclipid")], by=bloomevo[c("Cage", "phenotreatment","timepoint")], FUN=mean, na.rm=TRUE)
bloomevoap=bloomevomean[(bloomevomean$phenotreatment =="Apple"),]
bloomevobl=bloomevomean[(bloomevomean$phenotreatment =="Bloom"),]

appleevo=Evo[(Evo$cagetreatment =="Apple"),]
appleevomean <-aggregate(appleevo[c("dryweight","lipidweight", "arsqperclipid")], by=appleevo[c("Cage", "phenotreatment","timepoint")], FUN=mean, na.rm=TRUE)
appleevoap=appleevomean[(appleevomean$phenotreatment =="Apple"),]
appleevobl=appleevomean[(appleevomean$phenotreatment =="Bloom"),]

ATevo=Evo[(Evo$cagetreatment =="AT"),]
ATevomean <-aggregate(ATevo[c("dryweight","lipidweight", "arsqperclipid")], by=ATevo[c("Cage", "phenotreatment","timepoint")], FUN=mean, na.rm=TRUE)
ATevon=ATevomean[(ATevomean$phenotreatment =="Normal"),]
ATevox=ATevomean[(ATevomean$phenotreatment =="Axenic"),]

LBevo=Evo[(Evo$cagetreatment =="LB"),]
LBevomean <-aggregate(LBevo[c("dryweight","lipidweight", "arsqperclipid")], by=LBevo[c("Cage", "phenotreatment","timepoint")], FUN=mean, na.rm=TRUE)
LBevon=LBevomean[(LBevomean$phenotreatment =="Normal"),]
LBevox=LBevomean[(LBevomean$phenotreatment =="Axenic"),]




evomean <-aggregate(Evo[c("dryweight","lipidweight", "arsqperclipid")], by=Evo[c("cagetreatment", "phenotreatment","timepoint")], FUN=mean, na.rm=TRUE)
#bloomevomean=evomean[(evomean$`Cage Treatment` =="Bloom"),]
#appleevomean=evomean[(evomean$'Cage Treatment' =="Apple"),]
#LBevomean=evomean[(evomean$'Cage Treatment' =="LB"),]
#ATevomean=evomean[(evomean$'Cage Treatment'=="AT"),]
#APATLBevomeans=evomean[!(evomean$'Cage Treatment' == "Bloom"),]
evosummarybs <- ddply(Evo, c("cagetreatment","phenotreatment","timepoint"), summarise,
              N    = sum(!is.na(dryweight)),
              mean = mean(dryweight, na.rm=TRUE),
              sd   = sd(dryweight, na.rm=TRUE),
              se   = sd / sqrt(N)
)

#View(evosummarybs)

evosummarylipid <- ddply(Evo, c("cagetreatment","phenotreatment","timepoint"), summarise,
                      N    = sum(!is.na(lipidweight)),
                      mean = mean(lipidweight , na.rm=TRUE),
                      sd   = sd(lipidweight , na.rm=TRUE),
                      se   = sd / sqrt(N)
)

#View(evosummarylipid)

evosummaryperclipid <- ddply(Evo, c("cagetreatment","phenotreatment","timepoint"), summarise,
                         N    = sum(!is.na(arsqperclipid)),
                         mean = mean(arsqperclipid , na.rm=TRUE),
                         sd   = sd(arsqperclipid , na.rm=TRUE),
                         se   = sd / sqrt(N)
)

View(evosummaryperclipid)

####Subset Aggregate#### 
##lipid weight
bloomlipsum=subset(evosummarylipid,evosummarylipid$cagetreatment=="Bloom")
Applelipsum=subset(evosummarylipid,evosummarylipid$cagetreatment=="Apple")
ATlipsum=subset(evosummarylipid,evosummarylipid$cagetreatment=="AT")
LBlipsum=subset(evosummarylipid,evosummarylipid$cagetreatment=="LB")

APATLBlipsum=evosummarylipid[!(evosummarylipid$cagetreatment == "Bloom"),]
APATLBlipsumn=APATLBlipsum[!(APATLBlipsum$phenotreatment =="Axenic"),]
APATLBlipsumnp=APATLBlipsumn[!(APATLBlipsumn$phenotreatment =="Bloom"),]

AllEvosumlip=evosummarylipid[!(evosummarylipid$phenotreatment =="Axenic"),]
AllEvosumlipn=AllEvosumlip[!(AllEvosumlip$cagetreatment =="Apple" & AllEvosumlip$phenotreatment =="Bloom"),]
AllEvosumlipnp=AllEvosumlipn[!(AllEvosumlipn$cagetreatment =="Bloom" & AllEvosumlipn$phenotreatment =="Apple"),]


##dry weight 
bloombssum=subset(evosummarybs,evosummarybs$cagetreatment=="Bloom")
Applebssum=subset(evosummarybs,evosummarybs$cagetreatment=="Apple")
ATbssum=subset(evosummarybs,evosummarybs$cagetreatment=="AT")
LBbssum=subset(evosummarybs,evosummarybs$cagetreatment=="LB")


APATLBbssum=evosummarybs[!(evosummarybs$cagetreatment == "Bloom"),]
APATLBbssumn=APATLBbssum[!(APATLBbssum$phenotreatment =="Axenic"),]
APATLBbssumnp=APATLBbssumn[!(APATLBbssumn$phenotreatment =="Bloom"),]

AllEvosumbs=evosummarybs[!(evosummarybs$phenotreatment =="Axenic"),]
AllEvosumbsn=AllEvosumbs[!(AllEvosumbs$cagetreatment =="Apple" & AllEvosumbs$phenotreatment =="Bloom"),]
AllEvosumbsnp=AllEvosumbsn[!(AllEvosumbsn$cagetreatment =="Bloom" & AllEvosumbsn$phenotreatment =="Apple"),]

## percent lipid weight
bloomperclipsum=subset(evosummaryperclipid,evosummaryperclipid$cagetreatment=="Bloom")
Appleperclipsum=subset(evosummaryperclipid,evosummaryperclipid$cagetreatment=="Apple")
ATperclipsum=subset(evosummaryperclipid,evosummaryperclipid$cagetreatment=="AT")
LBperclipsum=subset(evosummaryperclipid,evosummaryperclipid$cagetreatment=="LB")

APATLBperclipsum=evosummaryperclipid[!(evosummaryperclipid$cagetreatment == "Bloom"),]
APATLBperclipsumn=APATLBperclipsum[!(APATLBperclipsum$phenotreatment =="Axenic"),]
APATLBperclipsumnp=APATLBperclipsumn[!(APATLBperclipsumn$phenotreatment =="Bloom"),]

AllEvosumperclip=evosummaryperclipid[!(evosummaryperclipid$phenotreatment =="Axenic"),]
AllEvosumperclipn=AllEvosumperclip[!(AllEvosumperclip$cagetreatment =="Apple" & AllEvosumperclip$phenotreatment =="Bloom"),]
AllEvosumperclipnp=AllEvosumperclipn[!(AllEvosumperclipn$cagetreatment =="Bloom" & AllEvosumperclipn$phenotreatment =="Apple"),]

####GRAPHS####

##body size graphs 

####bloomington dw####
bloombs <- ggplot(bloombssum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Dry Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple", "Bloomington"),
                      values = c("Red1", "Dark Blue")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  "4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
bloombs
####here####

this
bloombs2<-bloombs+geom_line(linetype= "dashed" , data=bloomevoap, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
bloombs3<-bloombs2+geom_line(linetype= "solid", data=bloomevobl, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
bloombs4=bloombs3+geom_point(data = Applebloomfoundermean, aes(x=timepoint, y=dryweight, group=cagetreatment, color=cagetreatment), size=4)
bloombs4
####apple dw####
applebs <- ggplot(Applebssum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Dry Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple", "Bloomington"),
                      values = c("Red1", "Dark Blue")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
                            #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
applebs

applebs2<-applebs+geom_line(linetype= "solid" , data=appleevoap, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
applebs2
applebs3<-applebs2+geom_line(linetype= "dashed", data=appleevobl, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
applebs4<-applebs3+geom_point(data = Applebloomfoundermean, aes(x=timepoint, y=dryweight, group=cagetreatment, color=cagetreatment), size=4)
applebs4

####AT dw####
ATbs <- ggplot(ATbssum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Dry Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Treatment", 
                      breaks = c("Normal","Axenic"),
                      labels = c("Normal AT", "Axenic AT"),
                      values = c("Purple", "Pink")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
ATbs
ATbs2<-ATbs+geom_line(linetype= "solid" , data=ATevon, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
ATbs3<-ATbs2+geom_line(linetype= "dashed", data=ATevox, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
ATbs4<-ATbs+geom_point(data = ATaxfoundermean, aes(x=timepoint, y=dryweight, group=cagetreatment, color=phenotreatment), size=4)
View(ATaxfoundermean)
ATbs4
####LB dw####
LBbs <- ggplot(LBbssum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Dry Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Treatment", 
                      breaks = c("Normal","Axenic"),
                      labels = c("Normal LB", "Axenic LB"),
                      values = c("Dark Green", "Light Green")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
LBbs
LBbs2<-LBbs+geom_line(linetype= "solid" , data=LBevon, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
LBbs3<-LBbs2+geom_line(linetype= "dashed", data=LBevox, aes(x=timepoint, y=dryweight, group=Cage), color="grey70", size=.5, alpha=.4)
LBbs4<-LBbs+geom_point(data = LBaxfoundermean, aes(x=timepoint, y=dryweight, group=cagetreatment, color=phenotreatment), size=4)
#View(LBaxfoundermean)
LBbs4

####GRAPH-ATLBAPPLE dw-024####
APATLBbssumnpn5=APATLBbssumnp[!(APATLBbssumnp$timepoint =="TP5"),]
APATLBbs1 <- ggplot(APATLBbssumnpn5, aes(x=timepoint, y= mean, group=cagetreatment))+
  geom_point(aes( color=cagetreatment),size=4)+
  geom_line(aes(color=cagetreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=cagetreatment), size=1)+
  ylab("Pooled Dry Weight")+
  xlab("Sample Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_x_discrete(labels=c("Initiation", "Summer", "Fall"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations"),
                      values = c("#848FA2","#058ED9","#CC2D35")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

APATLBbs1
APATLBbs2<-APATLBbs1+geom_point(data = AppleATLBfoundermean, aes(x=timepoint, y=dryweight, group=cagetreatment, color=cagetreatment), size=4)
APATLBbs3<-APATLBbs2+geom_line(linetype= "solid" , data=ATevon, aes(x=timepoint, y=dryweight, group=Cage), color="#CC2D35", size=.2, alpha=0.4)
APATLBbs4<-APATLBbs3+geom_line(linetype= "solid", data=LBevon, aes(x=timepoint, y=dryweight, group=Cage), color="#058ED9", size=.4, alpha=0.4)
APATLBbs5<-APATLBbs4+geom_line(linetype= "solid", data=appleevoap, aes(x=timepoint, y=dryweight, group=Cage), color="#848FA2", size=.2, alpha=0.4)
APATLBbs2
####everything dw####

ALLbs1 <- ggplot(AllEvosumbsnp, aes(x=timepoint, y= mean, group=cagetreatment))+
  geom_point(aes( color=cagetreatment),size=4)+
  geom_line(aes(color=cagetreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=cagetreatment), size=1)+
  ylab("Pooled Dry Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Bloom","Apple","LB", "AT"),
                      labels = c("Bloomington Cages","Apple Cages","LB Cages", "AT Cages"),
                      values = c("Dark Blue","Red1","Dark Green", "Purple")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

ALLbs1
ALLbs2<-ALLbs1+geom_point(data = AppleBloomATLBfoundermean, aes(x=timepoint, y=dryweight, group=cagetreatment, color=cagetreatment), size=4)
ALLbs2


##lipid weight graphs 

####bloomington lipid####
bloomlip <- ggplot(bloomlipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Lipid Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple", "Bloomington"),
                      values = c("Red1", "Dark Blue")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
                            "4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
bloomlip

bloomlip2<-bloomlip+geom_line(linetype= "dashed" , data=bloomevoap, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
bloomlip3<-bloomlip2+geom_line(linetype= "solid", data=bloomevobl, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
bloomlip4<-bloomlip3+geom_point(data = Applebloomfoundermean, aes(x=timepoint, y=lipidweight, group=cagetreatment, color=cagetreatment), size=4)

bloomlip4


####apple lipid####
applelip <- ggplot(Applelipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Lipid Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple", "Bloomington"),
                      values = c("Red1", "Dark Blue")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
applelip
applelip2<-applelip+geom_line(linetype= "solid" , data=appleevoap, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
applelip3<-applelip2+geom_line(linetype= "dashed", data=appleevobl, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
applelip4<-applelip3+geom_point(data = Applebloomfoundermean, aes(x=timepoint, y=lipidweight, group=cagetreatment, color=cagetreatment), size=4)

applelip4

####AT lipid####
ATlip <- ggplot(ATlipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Lipid Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Normal","Axenic"),
                      labels = c("Normal AT", "Axenic AT"),
                      values = c("Purple", "Pink")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
ATlip
ATlip2<-ATlip+geom_line(linetype= "solid" , data=ATevon, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
ATlip3<-ATlip2+geom_line(linetype= "dashed", data=ATevox, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
ATlip4<-ATlip3+geom_point(data = ATaxfoundermean, aes(x=timepoint, y=lipidweight, group=cagetreatment, color=phenotreatment), size=4)
ATlip4

####LB lipid####
LBlip <- ggplot(LBlipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Lipid Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Normal","Axenic"),
                      labels = c("Normal LB", "Axenic LB"),
                      values = c("Dark Green", "Light Green")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
LBlip
LBlip2<-LBlip+geom_line(linetype= "solid" , data=LBevon, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
LBlip3<-LBlip2+geom_line(linetype= "dashed", data=LBevox, aes(x=timepoint, y=lipidweight, group=Cage), color="grey70", size=.5, alpha=.4)
LBlip4<-LBlip3+geom_point(data = LBaxfoundermean, aes(x=timepoint, y=lipidweight, group=cagetreatment, color=phenotreatment), size=4)
LBlip4
View(ATaxfoundermean)

####ATLBAPPLE lipid##
APATLBlip1 <- ggplot(APATLBlipsumnp, aes(x=timepoint, y= mean, group=cagetreatment))+
  geom_point(aes( color=cagetreatment),size=4)+
  geom_line(aes(color=cagetreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=cagetreatment), size=1)+
  ylab("Pooled Lipid Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT"),
                      labels = c("Apple Cages","LB Cages", "AT Cages"),
                      values = c("Red1","Dark Green", "Purple")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)


APATLBlip2<-APATLBlip1+geom_point(data = AppleATLBfoundermean, aes(x=timepoint, y=lipidweight, group=cagetreatment, color=cagetreatment), size=4)
APATLBlip2

#APATLB2<-APATLB1+geom_line(linetype= "solid" , data=ATpign, aes(x=Time.Point, y=AVG, group=Cage.ID), color="Purple", size=.3, alpha=0.5)

#APATLB3<-APATLB2+geom_line(linetype= "solid", data=LBpign, aes(x=Time.Point, y=AVG, group=Cage.ID), color="Dark Green", size=.3, alpha=0.5)

#APATLB4<-APATLB3+geom_line(linetype= "solid", data=applepigap, aes(x=Time.Point, y=AVG, group=Cage.ID), color="Red1", size=.3, alpha=0.5)

####everything lipid####

ALLlip1 <- ggplot(AllEvosumlipnp, aes(x=timepoint, y= mean, group=cagetreatment))+
  geom_point(aes( color=cagetreatment),size=4)+
  geom_line(aes(color=cagetreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=cagetreatment), size=1)+
  ylab("Pooled Lipid Weight")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Bloom","Apple","LB", "AT"),
                      labels = c("Bloomington Cages","Apple Cages","LB Cages", "AT Cages"),
                      values = c("Dark Blue","Red1","Dark Green", "Purple")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

ALLlip2<-ALLlip1+geom_point(data = AppleBloomATLBfoundermean, aes(x=timepoint, y=lipidweight, group=cagetreatment, color=cagetreatment), size=4)

ALLlip2

####Percent Lipid####


####bloomington %lipid####
View(bloomperclipsum)
bloomperclip <- ggplot(bloomperclipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Percent Lipid")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple", "Bloomington"),
                      values = c("Red1", "Dark Blue")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
                            "4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
bloomperclip

bloomperclip2<-bloomperclip+geom_line(linetype= "dashed" , data=bloomevoap, aes(x=timepoint, y=perclipid, group=Cage), color="grey70", size=.5, alpha=.4)
bloomperclip3<-bloomperclip2+geom_line(linetype= "solid", data=bloomevobl, aes(x=timepoint, y=perclipid, group=Cage), color="grey70", size=.5, alpha=.4)
bloomperclip4<-bloomperclip3+geom_point(data = Applebloomfoundermean, aes(x=timepoint, y=perclipid, group=cagetreatment, color=cagetreatment), size=4)

bloomperclip4


####apple %lipid####
appleperclip <- ggplot(Appleperclipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Percent Lipid")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple", "Bloomington"),
                      values = c("Red1", "Dark Blue")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
appleperclip
appleperclip2<-appleperclip+geom_line(linetype= "solid" , data=appleevoap, aes(x=timepoint, y=perclipid, group=Cage), color="grey70", size=.5, alpha=.4)
appleperclip3<-appleperclip2+geom_line(linetype= "dashed", data=appleevobl, aes(x=timepoint, y=perclipid, group=Cage), color="grey70", size=.5, alpha=.4)
appleperclip4<-appleperclip3+geom_point(data = Applebloomfoundermean, aes(x=timepoint, y=perclipid, group=cagetreatment, color=cagetreatment), size=4)

appleperclip4

####AT lipid####
ATperclip <- ggplot(ATperclipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Percent Lipid")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Normal","Axenic"),
                      labels = c("Normal AT", "Axenic AT"),
                      values = c("Purple", "Pink")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
ATperclip
ATperclip2<-ATperclip+geom_line(linetype= "solid" , data=ATevon, aes(x=timepoint, y=arsqperclipid, group=Cage), color="grey70", size=.5, alpha=.4)
View(ATevon)
ATperclip3<-ATperclip2+geom_line(linetype= "dashed", data=ATevox, aes(x=timepoint, y=perclipid, group=Cage), color="grey70", size=.5, alpha=.4)
ATperclip4<-ATperclip+geom_point(data = ATaxfoundermean, aes(x=timepoint, y=arsqperclipid, group=cagetreatment, color=phenotreatment), size=4)
ATperclip4

####LB lipid####
LBperclip <- ggplot(LBperclipsum, aes(x=timepoint, y= mean, group=phenotreatment))+
  geom_point(aes( color=phenotreatment ),size=4)+
  geom_line(aes(color=phenotreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=phenotreatment), size=1)+
  ylab("Pooled Percent Lipid")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Phenotyping Diet", 
                      breaks = c("Normal","Axenic"),
                      labels = c("Normal LB", "Axenic LB"),
                      values = c("Dark Green", "Light Green")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)
LBperclip
LBperclip2<-LBperclip+geom_line(linetype= "solid" , data=LBevon, aes(x=timepoint, y=perclipid, group=Cage), color="grey70", size=.5, alpha=.4)
LBperclip3<-LBperclip2+geom_line(linetype= "dashed", data=LBevox, aes(x=timepoint, y=perclipid, group=Cage), color="grey70", size=.5, alpha=.4)
LBperclip4<-LBperclip+geom_point(data = LBaxfoundermean, aes(x=timepoint, y=arsqperclipid, group=cagetreatment, color=phenotreatment), size=4)
LBperclip4
View(ATaxfoundermean)

####GRAPH-ATLBAPPLE arsqperclipid -024####
APATLBperclip1 <- ggplot(APATLBperclipsumnp, aes(x=timepoint, y= mean, group=cagetreatment))+
  geom_point(aes( color=cagetreatment),size=4)+
  geom_line(aes(color=cagetreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=cagetreatment), size=1)+
  ylab("Pooled Percent Lipid (Arcsine Sqrt transformed)")+
  xlab("Sample Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_x_discrete(labels=c("Initiation", "Summer", "Fall"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","LB", "AT"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations"),
                      values = c("#848FA2","#058ED9","#CC2D35")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

APATLBperclip1
APATLBperclip2<-APATLBperclip1+geom_point(data = AppleATLBfoundermean, aes(x=timepoint, y=arsqperclipid, group=cagetreatment, color=cagetreatment), size=4)
APATLBperclip3<-APATLBperclip2+geom_line(linetype= "solid" , data=ATevon, aes(x=timepoint, y=arsqperclipid, group=Cage), color="#CC2D35", size=.3, alpha=0.5)
APATLBperclip4<-APATLBperclip3+geom_line(linetype= "solid", data=LBevon, aes(x=timepoint, y=arsqperclipid, group=Cage), color="#058ED9", size=.3, alpha=0.5)
APATLBperclip5<-APATLBperclip4+geom_line(linetype= "solid", data=appleevoap, aes(x=timepoint, y=arsqperclipid, group=Cage), color="#000000", size=.3, alpha=0.5)
APATLBperclip2

####everything lipid####

ALLperclip1 <- ggplot(AllEvosumperclipnp, aes(x=timepoint, y= mean, group=cagetreatment))+
  geom_point(aes( color=cagetreatment),size=4)+
  geom_line(aes(color=cagetreatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=cagetreatment), size=1)+
  ylab("Pooled Percent Lipid")+
  xlab("Time Point")+
  #scale_color_manual(values = c("Dark Blue", "Red"))+
  #scale_color_discrete( name = "Phenotyping Diet", labels = c("Apple", "Bloomington"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Bloom","Apple","LB", "AT"),
                      labels = c("Bloomington Cages","Apple Cages","LB Cages", "AT Cages"),
                      values = c("Dark Blue","Red1","Dark Green", "Purple")) +
  #scale_x_discrete(breaks=c("0","2","4","5"),
  #labels=c("Founder","Summer", "Fall", "Late Fall"))
  #scale_x_discrete(labels=c("0" = "Founder", "2" = "Summer",
  #"4" = "Fall", "5" = "Late Fall" ))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

ALLperclip2<-ALLperclip1+geom_point(data = AppleBloomATLBfoundermean, aes(x=timepoint, y=perclipid, group=cagetreatment, color=cagetreatment), size=4)

ALLperclip2

