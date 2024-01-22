##importing and librarying
library(ggplot2)
all= Orchard_2020_Phenotype_Master_528
library(tidyverse)
library(dplyr)
library(plyr)
#library(ddply)
library(cowplot)
library(nlme)
##Subsetting 
all.L=subset(all, all$Generation %in% c("F2", "F3"))

all.c=subset(all.L,all.L$Cage.Treatment %in% c("Start","A","B"))

View (all.c)
all.ma=subset(all,all$Cage.Treatment %in% c("start", "LB", "AT"))
all.L.sac=subset(all.L, all.L$Pheno.Treatment== "SAC")
all.L.0245=subset(all.L, all.L$Timepoint %in% c(0,2,4,5))
cage.pheno.f= factor(all.L.0245$Cage.Pheno)
all.L.sac.0245=subset(all.L.sac, all.L.sac$Timepoint %in% c(0,2,4,5))

###ANALYSIS LINEAR MIXED MODEL

all.L$cagextreatment=paste(all.L$Cage, all.L$Cage.Treatment)
all.c$cagextreatment=paste(all.c$Cage,all.c$Cage.Treatment)

## make summary data frame with only APPLE, BLOOM SAC
All.c.sum=all.c %>%
  group_by(Timepoint,Cage.Treatment,Cage) %>%
  #filter(!(Cage.Treatment == "B"))%>%
  filter(!any(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  filter(!(Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  #filter(!(Timepoint== 1))%>% ## remove timepoint 1
  #filter(!(Timepoint== 3))%>% ## remove timepoint 3
  #filter(!(Timepoint== 5))%>% ##remove timepoint 5
  mutate(cagetreatment = case_when(
    Cage.Treatment == "A" & Pheno.Treatment == "SAC" ~ "AA",
    Cage.Treatment == "A" & Pheno.Treatment == "B" ~ "AB",
    Cage.Treatment == "B" & Pheno.Treatment == "SAC" ~ "BB",
    Cage.Treatment == "B" & Pheno.Treatment == "A" ~ "BA", TRUE ~ "other"))

View(All.c.sum)

all.c.sum=All.c.sum %>%
  group_by(Timepoint,Cage.Treatment,Cage) %>%
  dplyr::summarize(N=n(),
            Cage=as.factor(Cage),
            cagetreatment=as.factor(cagextreatment),
            Cage.Treatment=as.factor(Cage.Treatment),
            Pheno.Treatment=as.factor(Pheno.Treatment),
            treatment=as.factor(Cage.Treatment),
            LD=mean(LD.Mean,na.rm = TRUE),
            LDm=mean(MLD.Mean,na.rm = TRUE),
            LDf=mean(FLD.Mean,na.rm = TRUE),
            V=mean(Viability,na.rm = TRUE),
            Starv=mean(Starv.Mean,na.rm = TRUE),
            Starvm=mean(MStarv.Mean,na.rm = TRUE),
            Starvf=mean(FStarv.Mean,na.rm =TRUE),
            Fecun=mean(F.1DayAVG,na.rm = TRUE),
  )
apbloomdata=unique(all.c.sum)

view(apbloomdata)
########lmm analysis with this subsetted dataframe

ABLmmdata=subset(apbloomdata, apbloomdata$Timepoint!='0') ##remove founder data
View(ABLmmdata)

LD1<-lme(LD ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(LD1) ##all SIG

LDm1<-lme(LDm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(LDm1) ##all SIG

LDf1<-lme(LDf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(LDf1) ##all SIG

V1<-lme(V ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(V1)   ##timepoint and intercept SIG , no treatment alone

Starv1<-lme(Starv ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(Starv1)  ##All SIG

Starvm1<-lme(Starvm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(Starvm1) ##All SIG

Starvf1<-lme(Starvf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(Starvf1) ##All SIG

Fecun1<-lme(Fecun ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=ABLmmdata)
anova(Fecun1) ##Treatment and Intercept SIG no Timepoint 


## make summary data frame with only APPLE and BLOOM common gardened
AB.c.sum=all.c %>%
  group_by(Timepoint,Cage.Treatment,Cage) %>%
  #filter(!(Cage.Treatment == "B"))%>%
  #filter(!any(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  #filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  #filter(!(Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  filter(!(Timepoint== 0))%>% 
  filter(!(Timepoint== 1))%>% ## remove timepoint 1
  filter(!(Timepoint== 3))%>% ## remove timepoint 3
  filter(!(Timepoint== 5))%>% ##remove timepoint 5
  mutate(cagetreatment = case_when(
    Cage.Treatment == "A" & Pheno.Treatment == "SAC" ~ "AA",
    Cage.Treatment == "A" & Pheno.Treatment == "B" ~ "AB",
    Cage.Treatment == "B" & Pheno.Treatment == "SAC" ~ "BB",
    Cage.Treatment == "B" & Pheno.Treatment == "A" ~ "BA", TRUE ~ "other"))

View(AB.c.sum)

AB.c.sum=AB.c.sum %>%
  group_by(Timepoint,cagetreatment,Cage) %>%
  dplyr::summarize(N=n(),
                   Cage=as.factor(Cage),
                   cagetreatment=as.factor(cagetreatment),
                   Cage.Treatment=as.factor(Cage.Treatment),
                   Pheno.Treatment=as.factor(Pheno.Treatment),
                   treatment=as.factor(cagetreatment),
                   LD=mean(LD.Mean,na.rm = TRUE),
                   LDm=mean(MLD.Mean,na.rm = TRUE),
                   LDf=mean(FLD.Mean,na.rm = TRUE),
                   V=mean(Viability,na.rm = TRUE),
                   Starv=mean(Starv.Mean,na.rm = TRUE),
                   Starvm=mean(MStarv.Mean,na.rm = TRUE),
                   Starvf=mean(FStarv.Mean,na.rm =TRUE),
                   Fecun=mean(F.1DayAVG,na.rm = TRUE),
  )
ABapbloomdata=unique(AB.c.sum)

Bapbloomdata=ABapbloomdata %>%
  filter(!(cagetreatment== "AA"))%>% 
  filter(!(cagetreatment== "BA"))

Aapbloomdata=ABapbloomdata %>%
  filter(!(cagetreatment== "BB"))%>% 
  filter(!(cagetreatment== "AB"))

view(Aapbloomdata)
########lmm analysis with this subsetted dataframe


BLD1<-lme(LD ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BLD1) 

ALD1<-lme(LD ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(ALD1) 


BLDm1<-lme(LDm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BLDm1)

ALDm1<-lme(LDm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(ALDm1) 


BLDf1<-lme(LDf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BLDf1) 

ALDf1<-lme(LDf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(ALDf1) 


BV1<-lme(V ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BV1)  

AV1<-lme(V ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(AV1) 


BStarv1<-lme(Starv ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BStarv1) 

AStarv1<-lme(Starv ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(AStarv1) 


BStarvm1<-lme(Starvm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BStarvm1) 

AStarvm1<-lme(Starvm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(AStarvm1) 


BStarvf1<-lme(Starvf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BStarvf1) 

AStarvf1<-lme(Starvf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(AStarvf1) 


BFecun1<-lme(Fecun ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Bapbloomdata)
anova(BFecun1) 

AFecun1<-lme(Fecun ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Aapbloomdata)
anova(AFecun1) 

#####################################
###FIGURE GENERATION
## regroup the same subset and make SEs
library(cowplot)
View(all.c)
all.c.sum=all.c %>%
  group_by(Timepoint, Cage.Treatment, Pheno.Treatment, Cage.Treatment.Pheno) %>%
  #filter(!(Cage.Treatment == "B"))%>%
  #filter(!(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  #filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  filter(!(Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  #filter(!(Timepoint== 1))%>% ## remove timepoint 1
  #filter(!(Timepoint== 3))%>% ## remove timepoint 3
  #filter(!(Timepoint== 5))%>% ##remove timepoint 5
  mutate(cagetreatment = case_when(
    Cage.Treatment == "A" & Pheno.Treatment == "SAC" ~ "AA",
    Cage.Treatment == "A" & Pheno.Treatment == "B" ~ "AB",
    Cage.Treatment == "B" & Pheno.Treatment == "SAC" ~ "BB",
    Cage.Treatment == "B" & Pheno.Treatment == "A" ~ "BA",
    TRUE                      ~ "other"))%>%
  mutate(Pheno.Treatment = case_when(
    Pheno.Treatment == "SAC" ~ "SAC",
    Pheno.Treatment == "A" ~ "Switched",
    Pheno.Treatment == "B" ~ "Switched"))%>%
  mutate(Cage.Treatment.Pheno = case_when(
    Cage.Treatment.Pheno == "BSAC" ~ "b",
    Cage.Treatment.Pheno == "ASAC" ~ "a",
    Cage.Treatment.Pheno == "AB" ~ "b",
    Cage.Treatment.Pheno == "BA" ~ "a",
    Cage.Treatment.Pheno == "StartSAC" ~ "b",
    Cage.Treatment.Pheno == "StartA" ~ "a"))%>%
  dplyr::summarize(N=n(),
            cagetreatment=as.factor(cagetreatment),
            phenotreatment=as.factor(Cage.Treatment.Pheno),
            pheno=as.factor(Pheno.Treatment),
            LD=mean(LD.Mean,na.rm = TRUE),
            LDsd=sd(LD.Mean,na.rm = TRUE),
            LDse= LDsd/sqrt(N),
            LDm=mean(MLD.Mean,na.rm = TRUE),
            LDmsd=sd(MLD.Mean,na.rm = TRUE),
            LDmse= LDmsd/sqrt(N),
            LDf=mean(FLD.Mean,na.rm = TRUE),
            LDfsd=sd(FLD.Mean,na.rm = TRUE),
            LDfse=LDfsd/sqrt(N),
            V=mean(Viability,na.rm = TRUE),
            V2=V*.85,
            Vsd=sd(Viability,na.rm = TRUE),
            Vse= Vsd/sqrt(N),
            Starv=mean(Starv.Mean,na.rm = TRUE),
            Starvsd=sd(Starv.Mean,na.rm = TRUE),
            Starvse= Starvsd/sqrt(N),
            Starvm=mean(MStarv.Mean,na.rm = TRUE),
            Starvmsd=sd(MStarv.Mean,na.rm = TRUE),
            Starvmse= Starvmsd/sqrt(N),
            Starvf=mean(FStarv.Mean,na.rm =TRUE),
            Starvfsd=sd(FStarv.Mean,na.rm = TRUE),
            Starvfse= Starvfsd/sqrt(N),
            Fecun=mean(F.1DayAVG,na.rm = TRUE),
            Fecunsd=sd(F.1DayAVG,na.rm = TRUE),
            Fecunse=Fecunsd/sqrt(N)
  )
apbloomfigdata=unique(all.c.sum)
View(apbloomfigdata)
apbloomfigdatanopheno=subset(apbloomfigdata, apbloomfigdata$Pheno.Treatment %in% c("SAC"))
View(all.c)


### for making individual cage lines
#allAPBLOOM<-subset(all.c,all.c$Pheno.Treatment=="SAC") 
allBloomonly<-subset(all.c, all.c$Cage.Treatment=="B" & all.c$Pheno.Treatment=="SAC") 
allAppleonly<-subset(all.c, all.c$Cage.Treatment=="A" & all.c$Pheno.Treatment=="SAC")

#LDapbloom<-aggregate(allAPBLOOM[c("LD.Mean")], by=allAPBLOOM[c("Timepoint","Cage","Cage.Treatment")], FUN=mean,)
aponly=aggregate(allAppleonly[c("LD.Mean", "MLD.Mean", "FLD.Mean", "Viability", "Starv.Mean", "FStarv.Mean", "MStarv.Mean", "F.1DayAVG" )], by=allAppleonly[c("Timepoint","Cage","Cage.Treatment")], FUN=mean,)
aponly$V2=(aponly$Viability*.85)

Bloomonly=aggregate(allBloomonly[c("LD.Mean", "MLD.Mean", "FLD.Mean", "Viability", "Starv.Mean", "FStarv.Mean", "MStarv.Mean", "F.1DayAVG")], by=allBloomonly[c("Timepoint","Cage","Cage.Treatment")], FUN=mean,)
Bloomonly$V2=(Bloomonly$Viability*.85)
View(aponly)

#SFgdata= subset(gdata, gdata$Timepoint != 0)

##LD Fig 
##no pheno treatment
LDg <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=LD, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment, shape=pheno),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=Cage.Treatment), size=1)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("A", "B", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()

LDg

##add individual cage lines

LDgl<-LDg+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=LD.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=LD.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))


LDgl

## with phenotype lines
#View(apbloomfigdata)
LDpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=LD, group=cagetreatment))+
  geom_point(aes(color=cagetreatment, shape=phenotreatment),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=cagetreatment), size=1)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_manual(name = "Cage Treatment", 
                      breaks = c("AA", "BB", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("red","black","red", "black")) +
  theme_cowplot()

LDpg


##test 
scale_custom= list(
  scale_color_manual(values = c("red","red","black","black"))
)


LDpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=LD, group=interaction(Cage.Treatment, phenotreatment)))+
  geom_point(aes( color=phenotreatment),size=4,show.legend = FALSE)+
  geom_line(aes(color=Cage.Treatment), show.legend = FALSE)+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=Cage.Treatment), size=.5, show.legend = FALSE)+
  geom_point(aes( color=Cage.Treatment),size=2, show.legend = FALSE)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_manual(values = c("red","red","black","black"))+
  theme_cowplot()

LDpg
##male LD fig ##
LDmg <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=LDm, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Cage.Treatment), size=1)+
  ylab("Male Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AA", "BB", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()

LDmg 

LDmgl<-LDmg+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=MLD.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=MLD.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))


LDmgl


LDmpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=LDm, group=cagetreatment))+
  geom_point(aes(color=cagetreatment, shape=pheno),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=cagetreatment), size=1)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AA", "BB", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()

LDmpg

## female LD figs##
LDfg <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=LDf, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Cage.Treatment), size=1)+
  ylab("Female Larval Development Time(hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("A", "B", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()

LDfg

LDfgl<-LDfg+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=FLD.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=FLD.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))


LDfgl

LDfpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=LDf, group=cagetreatment))+
  geom_point(aes(color=cagetreatment, shape=pheno),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=cagetreatment), size=1)+
  ylab("Female Larval Development Time(hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AA", "BB", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()

LDfpg

##viability figs##
Vg <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=V2, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=V2-Vse,ymax=V2+Vse, color=Cage.Treatment), size=1)+
  ylab("Egg to Adult Viability")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("A", "B", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()
Vg
Vgl<-Vg+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=V2, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=V2, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))
Vgl



Vpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=V2, group=interaction(Cage.Treatment, phenotreatment)))+
  geom_point(aes( color=phenotreatment),size=4,show.legend = FALSE)+
  geom_line(aes(color=Cage.Treatment), show.legend = FALSE)+
  geom_linerange(aes(ymin=V2-Vse,ymax=V2+Vse, color=Cage.Treatment), size=.5, show.legend = FALSE)+
  geom_point(aes( color=Cage.Treatment),size=2, show.legend = FALSE)+
  ylab("Viability")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_manual(values = c("red","red","black","black"))+
  theme_cowplot()

Vpg
## starvation figs##
STARVg <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=Starv, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=Starv-Starvse,ymax=Starv+Starvse, color=Cage.Treatment), size=1)+
  ylab("Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("A", "B", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()
STARVg 

STARVgl<-STARVg+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=Starv.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=Starv.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))
STARVgl



STARVpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=Starv, group=interaction(Cage.Treatment, phenotreatment)))+
  geom_point(aes( color=phenotreatment),size=4,show.legend = FALSE)+
  geom_line(aes(color=Cage.Treatment), show.legend = FALSE)+
  geom_linerange(aes(ymin=Starv-Starvse,ymax=Starv+Starvse, color=Cage.Treatment), size=.5, show.legend = FALSE)+
  geom_point(aes( color=Cage.Treatment),size=2, show.legend = FALSE)+
  ylab("Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_manual(values = c("red","red","black","black"))+
  theme_cowplot()

STARVpg

##male starvation figs##

STARVmg <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=Starvm, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=Starvm-Starvmse,ymax=Starvm+Starvmse, color=Cage.Treatment), size=1)+
  ylab(" Male Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("A", "B", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()
STARVmg 

STARVmgl<-STARVmg+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=MStarv.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=MStarv.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))
STARVmgl



STARVmg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=Starvm, group=cagetreatment))+
  geom_point(aes(color=cagetreatment, shape=pheno),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=Starvm-Starvmse,ymax=Starvm+Starvmse, color=cagetreatment), size=1)+
  ylab(" Male Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AA", "BB", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("red","black","red", "black")) +
  theme_cowplot()
STARVmg 

## female starv figs##

STARVfg <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=Starvf, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=Starvf-Starvfse,ymax=Starvf+Starvfse, color=Cage.Treatment))+
  ylab("Female Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("A", "B", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()
STARVfg 

STARVfgl<-STARVfg+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=FStarv.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=FStarv.Mean, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))
STARVfgl


STARVfpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=Starvf, group=cagetreatment))+
  geom_point(aes(color=cagetreatment, shape=pheno),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=Starvf-Starvfse,ymax=Starvf+Starvfse, color=cagetreatment))+
  ylab("Female Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("AA", "BB", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()
STARVfpg 

##fecund figs##
Fecung <- ggplot(apbloomfigdatanopheno, aes(x=as.factor(Timepoint), y=Fecun, group=Cage.Treatment))+
  geom_point(aes(color=Cage.Treatment),size=5)+
  geom_line(aes(color=Cage.Treatment))+
  geom_linerange(aes(ymin=Fecun-Fecunse,ymax=Fecun+Fecunse, color=Cage.Treatment), size=1)+
  ylab("Fecundity (Eggs/Day")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("A", "B", "AB", "BA"),
                      labels = c("Apple","Bloomington", "Apple on Bloomington", "Bloomington on Apple"),
                      values = c("#CC2D35","black","#CC2D35", "black")) +
  theme_cowplot()
Fecung 
View(aponly)
Fecungl<-Fecung+
  geom_line(data=aponly, aes(x=as.factor(Timepoint), y=F.1DayAVG, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  geom_line(data=Bloomonly, aes(x=as.factor(Timepoint), y=F.1DayAVG, group=as.factor(Cage), color=Cage.Treatment), 
            size=.5, alpha=.2) + 
  scale_color_manual(values=c(A="red", B="black"), labels=c("Apple", "Bloomington"))
Fecungl



FECUNpg <- ggplot(apbloomfigdata, aes(x=as.factor(Timepoint), y=Fecun, group=interaction(Cage.Treatment, phenotreatment)))+
  geom_point(aes( color=phenotreatment),size=4,show.legend = FALSE)+
  geom_line(aes(color=Cage.Treatment), show.legend = FALSE)+
  geom_linerange(aes(ymin=Fecun-Fecunse,ymax=Fecun+Fecunse, color=Cage.Treatment), size=.5, show.legend = FALSE)+
  geom_point(aes( color=Cage.Treatment),size=2, show.legend = FALSE)+
  ylab("Fecundity (Eggs / Day)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_manual(values = c("red","red","black","black"))+
  theme_cowplot()

FECUNpg
####NEED PIGMENTATION






#####################################################################################################################################

###adding individual cage lines 
View(apbloomfigdata)
##Bloomington only 
Bloomfigdata=subset(apbloomfigdata, apbloomfigdata$Cage.Treatment== "B" & apbloomfigdata$Pheno.Treatment== "SAC")
Applefigdata=subset(apbloomfigdata, apbloomfigdata$Cage.Treatment== "A" & apbloomfigdata$Pheno.Treatment== "SAC")
BloomAppplefigdata=subset(apbloomfigdata, apbloomfigdata$Pheno.Treatment=="SAC")

View(Bloomfigdata)
BloomLD <- ggplot(Bloomfigdata, aes(x=as.factor(Timepoint), y=LD, group=cagetreatment))+
  geom_point(aes(color=cagetreatment),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=cagetreatment), size=1)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_discrete(name = "Cage Treatment") +
  theme_cowplot()
BloomLD 


allAPBLOOM<-subset(all.c,all.c$Pheno.Treatment=="SAC") 
allBloomonly<-subset(all.c, all.c$Cage.Treatment=="B" & all.c$Pheno.Treatment=="SAC") 
allAppleonly<-subset(all.c, all.c$Cage.Treatment=="A" & all.c$Pheno.Treatment=="SAC")

LDapbloom<-aggregate(allAPBLOOM[c("LD.Mean")], by=allAPBLOOM[c("Timepoint","Cage","Cage.Treatment")], FUN=mean, na.rm=TRUE)

BLOOMLD<-BloomLD+geom_line(data=LDapbloom, aes(x=Timepoint, y=LD.Mean, group=Cage, color=Cage.Treatment), size=.5)

BLOOMLD

AppleLD <- ggplot(Applefigdata, aes(x=as.factor(Timepoint), y=LD, group=cagetreatment))+
  geom_point(aes(color=cagetreatment),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=cagetreatment), size=1)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_discrete(name = "Cage Treatment") +
  theme_cowplot()

LDApple<-aggregate(allAppleonly[c("LD.Mean")], by=allAppleonly[c("Timepoint","Cage")], FUN=mean, na.rm=TRUE)
LDApple$Cage<-factor(LDApple$Cage, levels=c("B","A1", "A2", "A3", "A4", "A5", "A6"))
LDApple$Timepoint<-factor(LDApple$Timepoint, levels=c("0", "1", "2", "3", "4", "5"))

APPLELD<-AppleLD+geom_line(data=LDApple, aes(x=Timepoint, y=LD.Mean, group=Cage), color="grey70", size=.5)
APPLELD 

BothLD=ggplot(BloomAppplefigdata, aes(x=as.factor(Timepoint), y=LD, group=cagetreatment))+
  geom_point(aes(color=cagetreatment),size=5)+
  geom_line(aes(color=cagetreatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=cagetreatment), size=1)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "TP1", "TP2", "TP3", "TP4", "TP5"))+
  scale_color_discrete(name = "Cage Treatment") +
  theme_cowplot()

BOTHLD=BothLD+
  geom_line(data=LDbloom, aes(x=Timepoint, y=LD.Mean, group=Cage), color="grey70", size=.5)+
  geom_line(data=LDApple, aes(x=Timepoint, y=LD.Mean, group=Cage), color="grey70", size=.5)

BOTHLD
