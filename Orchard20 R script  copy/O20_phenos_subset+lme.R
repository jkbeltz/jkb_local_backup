##importing and librarying
library(ggplot2)
all= Orchard_2020_Phenotype_Master_528
library(tidyverse)
library(dplyr)
library(cowplot)
library(nlme)
library(lme4)
##Subsetting 
all.L=subset(all, all$Generation %in% c("F2", "F3"))
all.c=subset(all,all$Cage.Treatment %in% c("Start","A","B"))
all.ma=subset(all,all$Cage.Treatment %in% c("start", "LB", "AT"))
all.L.sac=subset(all.L, all.L$Pheno.Treatment== "SAC")
all.L.0245=subset(all.L, all.L$Timepoint %in% c(0,2,4,5))
cage.pheno.f= factor(all.L.0245$Cage.Pheno)
all.L.sac.0245=subset(all.L.sac, all.L.sac$Timepoint %in% c(0,2,4,5))

###ANALYSIS 

all.L$cagextreatment=paste(all.L$Cage, all.L$Cage.Treatment)

## make summary data frame with only APPLE, LB, AT, SAC 024
all.L.sum=all.L %>%
  group_by(Timepoint, cagextreatment) %>%
  #filter(!(Cage.Treatment == "B"))%>%
  filter(!any(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  filter(!(Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  filter(!(Timepoint== 1))%>% ## remove timepoint 1
  filter(!(Timepoint== 3))%>% ## remove timepoint 3
  filter(!(Timepoint== 5))%>% ##remove timepoint 5
  summarise(N=n(),
            Cage=as.factor(Cage),
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
lmmdata=unique(all.L.sum)
View(lmmdata)

########lmm analysis with this subsetted dataframe

Lmmdata=drop.levels(subset(lmmdata, Timepoint!='0')) ##remove founder data

LD1<-lme(LD ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(LD1) ##treatment SIG

LDm1<-lme(LDm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(LDm1) ##treatment SIG

LDf1<-lme(LDf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(LDf1) ##treatment SIG

V1<-lme(V ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(V1)   ##timepoint SIG

Starv1<-lme(Starv ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(Starv1)  ##timepoint SIG

Starvm1<-lme(Starvm ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(Starvm1) ##timepoint SIG, treatment SIG, interaction ALMOST SIG

Starvf1<-lme(Starvf ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(Starvf1) ##timepoint SIG, interaction SIG

Fecun1<-lme(Fecun ~ Timepoint*treatment, random=~1|Cage/Timepoint, data=Lmmdata)
anova(Fecun1) ##NO SIG
?lme
  