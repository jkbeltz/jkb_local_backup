####importing and librarying####
library(ggplot2)
all= Orchard_2020_Phenotype_Master_This
library(tidyverse)
library(dplyr)
####Subsetting ####
all.L=subset(all, all$Generation %in% c("F2", "F3"))
View(all.L)
all.c=subset(all,all$Cage.Treatment %in% c("Start","A","B"))
all.ma=subset(all,all$Cage.Treatment %in% c("start", "LB", "AT"))
all.L.sac=subset(all.L, all.L$Pheno.Treatment== "SAC")
all.L.0245=subset(all.L, all.L$Timepoint %in% c(0,2,4))
cage.pheno.f= factor(all.L.0245$Cage.Pheno)
all.L.sac.0245=subset(all.L.sac, all.L.sac$Timepoint %in% c(0,2,4,5))

####ANALYSIS ####

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
            timepoint=as.factor(Timepoint),
            cageID=as.factor(Cage),
            cagetreatment=as.factor(Cage.Treatment),
            phenotreatment=as.factor(Pheno.Treatment),
            LD=mean(LD.Mean,na.rm = TRUE),
            LDMale=mean(MLD.Mean,na.rm = TRUE),
            LDFemale=mean(FLD.Mean,na.rm = TRUE),
            Viability=mean(Viability,na.rm = TRUE),
            Starv=mean(Starv.Mean,na.rm = TRUE),
            StarvMale=mean(MStarv.Mean,na.rm = TRUE),
            StarvFemale=mean(FStarv.Mean,na.rm =TRUE),
            Fecun=mean(F.1DayAVG,na.rm = TRUE),
  )
lmmdata=unique(all.L.sum)
View(Lmmdata)

########lmm analysis with this subsetted dataframe####

Lmmdata=subset(lmmdata, Timepoint!='0') ##remove founder data

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


####FIGURE GENERATION####
## regroup the same subset and make SEs
View(all.L)
all.L.sum=all.L %>%
  group_by(Timepoint, Cage.Treatment) %>%
  filter(!(Cage.Treatment == "B"))%>%
  filter(!(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  filter(!(Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  filter(!(Timepoint== 1))%>% ## remove timepoint 1
  filter(!(Timepoint== 3))%>% ## remove timepoint 3
  filter(!(Timepoint== 5))%>% ##remove timepoint 5
  summarise(N=n(),
            treatment=as.factor(Cage.Treatment),
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
gdata=unique(all.L.sum)
View(gdata)

SFgdata= subset(gdata, gdata$Timepoint != 0)

LDg <- ggplot(gdata, aes(x=as.factor(Timepoint), y=LD, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=treatment), size=2)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot()+
  theme(legend.position = "none")

LDg + theme(legend.position = c(0.6, 0.15))

SFLDg <- ggplot(SFgdata, aes(x=as.factor(Timepoint), y=LD, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=treatment), size=2)+
  ylab("Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c( "Summer", "Fall"))+
  scale_color_discrete(name = NULL, labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot()+
  theme(legend.position = "none")

SFLDg 



LDmg <- ggplot(gdata, aes(x=as.factor(Timepoint), y=LDm, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=treatment), size=2)+
  ylab("Male Larval Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.8,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

LDmg + theme(legend.position = c(0.6, 0.15))

LDfg <- ggplot(gdata, aes(x=as.factor(Timepoint), y=LDf, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=treatment), size=2)+
  ylab("Female Larval Development Time(hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.8,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

LDfg + theme(legend.position = c(0.6, 0.15))

Vg <- ggplot(gdata, aes(x=as.factor(Timepoint), y=V, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=V-Vse,ymax=V+Vse, color=treatment), size=2)+
  ylab("Egg to Adult Viability")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.8,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

Vg + theme(legend.position = c(0.65, 0.2))

STARVg <- ggplot(gdata, aes(x=as.factor(Timepoint), y=Starv, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=Starv-Starvse,ymax=Starv+Starvse, color=treatment), size=2)+
  ylab("Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.8,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

STARVg + theme(legend.position = c(0.65, 0.8))



STARVmg <- ggplot(gdata, aes(x=as.factor(Timepoint), y=Starvm, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=Starvm-Starvmse,ymax=Starvm+Starvmse, color=treatment), size=2)+
  ylab(" Male Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.8,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

STARVmg + theme(legend.position = c(0.65, 0.8))

STARVfg <- ggplot(gdata, aes(x=as.factor(Timepoint), y=Starvf, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=Starvf-Starvfse,ymax=Starvf+Starvfse, color=treatment), size=2)+
  ylab("Female Starvation Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.8,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

STARVfg + theme(legend.position = c(0.65, 0.8))


Fecung <- ggplot(gdata, aes(x=as.factor(Timepoint), y=Fecun, group=treatment))+
  geom_point(aes(color=treatment),size=5)+
  geom_line(aes(color=treatment))+
  geom_linerange(aes(ymin=Fecun-Fecunse,ymax=Fecun+Fecunse, color=treatment), size=2)+
  ylab("Fecundity (Eggs/Day")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.8,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

Fecung + theme(legend.position = c(0.65, 0.8)) +scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))


####Including sterilized samples####

all.L.sum2=all.L.0245%>%
  group_by(Timepoint, Cage.Treatment, Pheno.Treatment) %>%
  filter(!(Cage.Treatment == "B"))%>%
  filter(!(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  #filter(!(Timepoint < 1 & Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  filter(!(Timepoint== 1))%>% ## remove timepoint 1
  filter(!(Timepoint== 3))%>% ## remove timepoint 3
  filter(!(Timepoint== 5))%>% ##remove timepoint 5
  summarise(N=n(),
            treatment=as.factor(Cage.Treatment),
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
g2data=unique(all.L.sum2)
g2data$treatmentpheno=paste(g2data$Cage.Treatment,g2data$Pheno.Treatment)
View(g2data)


####replotting with sterilized samples####
library(cowplot)
LD2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=LD, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=Cage.Treatment), size=1)+
  ylab("Larval Developement Time (hrs)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Initiation ", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

LD2 + theme_cowplot()

LDm2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=LDm, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Cage.Treatment), size=1)+
  ylab("Male Larval Developement Time (hrs)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

LDm2 + theme_cowplot()

LDf2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=LDf, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Cage.Treatment), size=1)+
  ylab("Female Larval Developement Time (hrs)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

LDf2 + theme_cowplot()

Starv2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=Starv, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=Starv-Starvse,ymax=Starv+Starvse, color=Cage.Treatment), size=1)+
  ylab("Dessication Time (hrs)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

Starv2 + theme_cowplot()

Starvm2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=Starvm, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=Starvm-Starvmse,ymax=Starvm+Starvmse, color=Cage.Treatment), size=1)+
  ylab("Male Dessication Time (hrs)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

Starvm2 + theme_cowplot()

Starvf2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=Starvf, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=Starvf-Starvfse,ymax=Starvf+Starvfse, color=Cage.Treatment), size=1)+
  ylab("Female Starvation Time (hrs)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Initiation", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

Starvf2 + theme_cowplot() + theme(legend.position = c(0.64, 0.82))+ theme(legend.text=element_text(size=9))

SFg2data= subset(g2data, g2data$Timepoint != 0)
View(SFg2data)

SFStarvf2 <- ggplot(SFg2data, aes(x=as.factor(Timepoint), y=Starvf, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=Starvf-Starvfse,ymax=Starvf+Starvfse, color=Cage.Treatment), size=1)+
  ylab("Female Starvation Time (hrs)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))+
  theme_cowplot()+
  theme(legend.position = "none")

SFStarvf2 + theme_cowplot() + theme(legend.position = c(0.64, 0.88))+ theme(legend.text=element_text(size=9))


V2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=V, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=V-Vse,ymax=V+Vse, color=Cage.Treatment), size=1)+
  ylab("Egg to Adult Viability (%)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

V2 + theme_cowplot()

Fecun2 <- ggplot(g2data, aes(x=as.factor(Timepoint), y=Fecun, group=treatmentpheno))+
  geom_point(aes(shape=Pheno.Treatment, color=Cage.Treatment),size=4)+
  geom_line(aes(color=Cage.Treatment, linetype=Pheno.Treatment))+
  geom_linerange(aes(ymin=Fecun-Fecunse,ymax=Fecun+Fecunse, color=Cage.Treatment), size=1)+
  ylab("Fecundity (Eggs/Fly/Day)")+
  xlab("Sample Timepoint")+
  scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  scale_linetype(name=NULL, labels = c("Microbiome Present", "Sterilized"))+
  scale_shape_discrete(name=NULL, labels = c("Microbiome Present", "Sterilized")) +
  scale_colour_discrete(name= "Cage Treatment", labels = c("Apple Food", "Apple + Ac. Thailandicus", "Apple + Lac. Brevis"))

Fecun2 + theme_cowplot()

####founder####
all.L.sumfo=all.L %>%
  group_by(Timepoint, Cage.Treatment, Pheno.Treatment) %>%
  filter(!(Timepoint > 1 )) %>% #Removing apple treatments from bloomington cages
  filter(!(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing bloomington treatments from apple cages
  filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  filter(!(Timepoint== 1))%>% ## remove timepoint 1
  summarise(N=n(),
            treatment=as.factor(Cage.Treatment),
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
fdata=unique(all.L.sumfo)
fdata$treatmentpheno=paste(fdata$Cage.Treatment,fdata$Pheno.Treatment)
view(fdata)

fdata$treatmentpheno= factor(c( "Apple", "Sterilized Apple", "Apple + Ac. Thailandicus","Bloomington", "Apple + Lac. Brevis"))
colorset= c('Apple'='#F8766D', 'Sterilized Apple'= 'dark grey', 'Apple + Ac. Thailandicus' = '#00BA38','Bloomington'= 'purple', 'Apple + Lac. Brevis' = '#619CFF')

FBARLD <- ggplot(fdata, aes(x=treatmentpheno, y=LD))+
  scale_colour_manual(name = NULL, values = colorset)+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=treatmentpheno), size=15)+
  geom_point(pch=18, size=4)+
  ylab("Founder Larval Development Time (Hrs)")+
  xlab("Phenotyping Treatment")+
  theme_cowplot()+ 
  theme(axis.text.x = element_blank())+
  theme(legend.position = "none")

FBARLD 

FBARSTARV <- ggplot(fdata, aes(x=treatmentpheno, y=Starv))+
  scale_colour_manual(name = NULL, values = colorset)+
  geom_linerange(aes(ymin=Starv-Starvse,ymax=Starv+Starvse, color=treatmentpheno), size=15)+
  geom_point(pch=18, size=4)+
  ylab("Founder Starvation Time (Hrs)")+
  xlab("Phenotyping Treatment")+
  theme_cowplot()+
  theme(axis.text.x = element_blank())+
  theme(legend.position = "none")

FBARSTARV 

FBARFECUN <- ggplot(fdata, aes(x=treatmentpheno, y=Fecun))+
  geom_point(size=4)+
  scale_colour_manual(name = NULL, values = colorset)+
  geom_linerange(aes(ymin=Fecun-Fecunse,ymax=Fecun+Fecunse, color=treatmentpheno), size=1)+
  ylab("Founder Fecundity")+
  xlab("Phenotyping Treatment")

FBARFECUN + theme_cowplot() + theme(axis.text.x = element_blank())

####foudner compound figure####
install.packages("ggpubr")
library("ggpubr")

 ggarrange(
  FBARLD, FBARSTARV, SFLDg, SFStarvf2
)

 Found=ggarrange(
   FBARLD, FBARSTARV, common.legend = TRUE, legend = "right", labels= c("A", "B"))
Found

prelim=ggarrange(SFLDg, SFStarvf2, labels = c("A", "B"),common.legend = TRUE, legend = "right")
prelim
 
p =  p + guides(shape = guide_legend(override.aes = list(size = 0.5)))
p
p <- p + guides(color = guide_legend(override.aes = list(size = 0.5)))
p <- p + theme(legend.title = element_text(size = 3), 
               legend.text = element_text(size = 3))
  
p 
  
  
Treat <-factor(all.L.sum2$treat)
Pheno <-factor(all.L.sum2$pheno)

p=ggplot(all.L.sum2) +
  geom_point(aes(x = Timepoint, y= LD, color=Treat, shape=Pheno), 
             size = 2, alpha = .75, position=position_jitter(h=0.0,w=0.0))

p
p1<-p + geom_path( aes(x = Timepoint, y= LD, group=Cage.Treatment.Pheno, color=Treat))
p1


#######
#L.0245.sum.less = L.0245.sum[L.0245.sum$Cage.Pheno != "A1B","A2B", "A3B", "A4B", 
#                              "A5B", "A6B", "B1A","B2A",
#                              "B3A","B4A","B5A","B6A", , drop=FALSE]
#View(L.0245.sum.less)

#library(RColorBrewer)
#myColors <- brewer.pal(6,"Set1")
#names(myColors) <- levels(L.0245.sum2$treat)
#colScale <- scale_colour_manual(name = "treat",values = myColors)


Treat <-factor(L.0245.sum2$treat)
Pheno <-factor(L.0245.sum2$pheno)

p=ggplot(L.0245.sum2) +
    geom_point(aes(x = Timepoint, y= LD, color=Treat, shape=Pheno), 
               size = 2, alpha = .75, position=position_jitter(h=0.0,w=0.2))

p

Treat <-factor(L.0245.sum2$treat)
Pheno <-factor(L.0245.sum2$pheno)

p1<-p + geom_path( aes(x = Timepoint, y= LD, group=cage, color=Treat))


p1

pc <- p + colScale
pc
p1 <- p %+% droplevels(subset(L.0245.sum, L.0245.sum$Cage.Pheno != "A1B","A2B", "A3B", "A4B", 
                                                           "A5B", "A6B", "B1A","B2A",
                                                           "B3A","B4A","B5A","B6A"))

p1
tab %>%
    group_by(month, variable) %>%
    summarise(a_sum=sum(amount),
              a_mean=(mean(amount)))
View(data)
 #%>%
  group_by(month, variable) %>%
  summarise(a_sum=sum(amount),
            a_mean=(mean(amount)))

##Data Viz

a=ggplot(data=all.L.sac.0245, aes(x =as.factor(Timepoint), y=Viability))+ 
    xlab('Time Point') + ylab("Starvation Time") + 
    geom_boxplot(aes(fill=as.factor(Cage.Treatment.Pheno)))

