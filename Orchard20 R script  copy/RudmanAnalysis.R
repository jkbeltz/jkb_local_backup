### Made by Seth Rudman for E cage phnenotyping data
#SMR_PhenotypeScript
#Load some packages - 'install.packages('xxx')' Not all of these are needed here, but they are all useful
#install.packages("vegan")
library(vegan)
library(ecodist)
library(reshape2)
library(ggplot2)
library(gdata)
library(plyr)
library(compute.es)
library(lme4)
library(nlme)
library(car)
#install.packages("MASS", "cowplot")
#install.packages("cowplot")
library(MASS)
library(ggplot2)
library(cowplot)
library(tidyr)
library(multcomp)

#Analysis and plotting for phenotypes that are single measures (e.g., fecundity, body size, etc.)

#import data - make sure it's a .csv
#fecund <-read.csv(file= "~/Documents/DrosEED14/DrosEED_fecundity_170227.csv", stringsAsFactors= TRUE, strip.white = TRUE, na.strings = c("NA","na"))
fecund = fecundR
str(fecund)
## 'data.frame':    1610 obs. of  7 variables:
##  $ treatment  : Factor w/ 2 levels "Non-Seasonal",..: 2 2 2 2 2 2 2 2 2 2 ...
##  $ cage       : Factor w/ 20 levels "E1","E10","E2",..: 1 3 4 5 6 7 8 9 10 2 ...
##  $ sample     : Factor w/ 5 levels "founder","four",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ fecund1    : int  NA NA NA NA NA NA NA NA NA NA ...
##  $ fecund2    : int  NA NA NA NA NA NA NA NA NA NA ...
##  $ fecund3    : int  NA NA NA NA NA NA NA NA NA NA ...
##  $ totalfecund: int  NA NA NA NA NA NA NA NA NA NA ...
tail(fecund)
##      treatment cage sample fecund1 fecund2 fecund3 totalfecund
## 1605  Seasonal  E10   four       2      33      64          99
## 1606  Seasonal  E10   four       3      60      72         135
## 1607  Seasonal  E10   four       4      22      56          82
## 1608  Seasonal  E10   four       1      29      68          98
## 1609  Seasonal  E10   four      22      77     102         201
## 1610  Seasonal  E10   four       0      17      92         109
#change timepoint to a factor
fecund$sample<-factor(fecund$sample, levels=c("founder","one", "two", "three", "four", "five"))

##MAKE COLUMNS NUMERIC
View(fecund)
sapply(fecund, class)
cols.num <- c("fecund1","fecund2", "fecund3", "totalfecund")
cols.factor <- c("treatment", "sample", "cage", "assayenv", "pheno")
fecund[cols.num] <- sapply(fecund[cols.num],as.numeric)
fecund[cols.factor] <- sapply(fecund[cols.factor],as.factor)
sapply(fecund, class)
## Make perday column
fecund$perday<-fecund$totalfecund/3

##subset data for first analysis only apple, LB, AT/ SAC/ 0,2,4
fecund=subset(fecund, fecund$pheno == "same")
fecund=subset(fecund, fecund$treatment %in% c("Apple", "LB", "AT"))
fecund=subset(fecund, fecund$sample %in% c("founder", "two", "four"))
View(fecund)

#ggplot graph
##aggregate means around cage and treatment
fecmean <-aggregate(fecund[c("perday")], by=fecund[c("cage", "treatment", "sample")], FUN=mean, na.rm=TRUE)
##summary / find SE
cdat <- ddply(fecmean, c("treatment",  "sample"), summarise,
              N    = sum(!is.na(perday)),
              mean = mean(perday, na.rm=TRUE),
              sd   = sd(perday, na.rm=TRUE),
              se   = sd / sqrt(N)
)

cdat$sample<-factor(cdat$sample, levels=c("founder", "two", "four"))

#cdat1<-drop.levels(subset(cdat, treatment!="Non-Seasonal"))
colors<-c('black',  'grey70')
##clarify sample levels 
cdat1$sample<-factor(cdat1$sample, levels=c("founder", "two", "four"))
##plot summary data with COWPLOT VISUALS
g1 <- ggplot(cdat1, aes(x=sample, y=mean, group=treatment))+
  geom_point(aes( color=treatment),size=4)+
  geom_line(aes(color=treatment), size=1.5)+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=treatment), size=1)+
  ylab("Fecundity (eggs/day)")+
  xlab("Sample Point")+
  scale_x_discrete(labels=c("Founder","Summer", "Fall"))+
  theme_cowplot(font_size = 14, font_family = "", line_size = 0.5,
                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14)

g1

#Add profile of each individual cage
#subset to just 'E' cages
fecundE<-droplevels(subset(fecund, treatment=="Seasonal"))

fecE <-aggregate(fecundE[c("perday")], by=fecundE[c("treatment","sample", "cage")], FUN=mean, na.rm=TRUE)
fecE$cage<-factor(fecE$cage, levels=c("founder","E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10"))

gFec2<-g1+geom_line(data=fecE, aes(x=sample, y=perday, group=cage), color="grey70", size=.5)

gFec2

#LME based analysis for effects of time and treatment with cage as a random effect

#
#Fecundity analysis
#
#
fecund <-read.csv(file= "~/Documents/DrosEED14/DrosEED_fecundity_170227.csv", stringsAsFactors= TRUE, strip.white = TRUE, na.strings = c("NA","na"))
head(fecund)
##   treatment cage  sample fecund1 fecund2 fecund3 totalfecund
## 1  Seasonal   E1 founder      NA      NA      NA          NA
## 2  Seasonal   E2 founder      NA      NA      NA          NA
## 3  Seasonal   E3 founder      NA      NA      NA          NA
## 4  Seasonal   E4 founder      NA      NA      NA          NA
## 5  Seasonal   E5 founder      NA      NA      NA          NA
## 6  Seasonal   E6 founder      NA      NA      NA          NA
fecund$sample<-factor(fecund$sample, levels=c("founder","one", "two", "three", "four"))
fecund$perday<-fecund$totalfecund/3

#ggplot graph
fecmean <-aggregate(fecund[c("perday")], by=fecund[c("cage", "treatment", "sample")], FUN=mean, na.rm=TRUE)
fecmean1<-drop.levels(subset(fecmean, sample!="founder"))

#LME
x1<-lme(perday ~ sample*treatment, random=~1|cage/sample, data=fecmean1)
#summary(x1)
anova(x1)

#Analysis and plotting for phenotypes that are single measures (e.g., fecundity, body size, etc.) This code is based on data collect without individual entries (e.g., 5 flies dead at Tn, 7 flies dead at Tn+1) Uses a GLM to estimate the time at which half the flies developed/starved/knocked down/etc.
## repeaty for viability data 

viabile=ViabileR
viabile=subset(viabile, viabile$pheno == "same")
viabile=subset(viabile, viabile$treatment %in% c("Apple", "LB", "AT"))
viabile=subset(viabile, viabile$sample %in% c("founder", "two", "four"))
viabile$sample<-factor(viabile$sample, levels=c("founder", "two", "four"))
viamean <-aggregate(viabile[c("viability")], by=viabile[c("cage", "treatment", "sample")], FUN=mean, na.rm=TRUE)
viamean1<-drop.levels(subset(viamean, sample!="founder"))
View(viamean1)
#LME
x1<-lme(viability ~ sample*treatment, random=~1|cage/sample, data=viamean1)
#summary(x1)
anova(x1)

#Function to calculate LD50 - time at which 50 percent of individuals in a vial died/eclosed/etc.
f<- function(x) {
  fit<-glm(cbind(alive, value) ~ time, family="binomial", data=x)
  dose.p(fit, p=0.5)[[1]]
}

#
#
#Profiles of trait evolution over time for EED14 using LD50
#

#Dessication - data from founder
des0 <-read.csv(file= "~/Documents/DrosEED14/Dessication_DrosEED14_T0_180416.csv", stringsAsFactors= TRUE, strip.white = TRUE, na.strings = c("NA","na"), check.names= FALSE)
head(des0)

desall= desR

desall=subset(desall, desall$pheno == "same")
desall=subset(desall, desall$treatment %in% c("Apple", "LB", "AT"))
desall=subset(desall, desall$sample %in% c("founder", "two", "four"))
desall$sample<-factor(desall$sample, levels=c("founder", "two", "four"))

des0=subset(desall, desall$sample == "founder")
drop <- c("tubesample","cagesex", "treatmentsex", "replicate", "sample", "cagesample", "assayenv", "pheno")
des0 = des0[,!(names(des0) %in% drop)]
head(des0)

##   treatment    cage vial sex 0 3 6 8 10 12 14 16 18 20 22 24 26
## 1   founder founder   1M   M 0 0 0 0  4  7  8  8  8  8 10 10 10
## 2   founder founder   1F   F 0 0 0 0  0  1  3  3  7  9 10 10 10
## 3   founder founder   2M   M 0 0 0 0  2  6  9  9 10 10 10 10 10
## 4   founder founder   2F   F 0 0 0 0  0  0  2  4  4  6 10 10 10
## 5   founder founder   3M   M 0 0 0 0  1  3  9 10 10 10 10 10 10
## 6   founder founder   3F   F 0 0 0 0  0  0  1  5  7  9 10 10 10
des0m<-melt(des0)
View(des0m)
## Using treatment, cage, vial, sex as id variables
names(des0m) = c("cage","treatment","sex","vial", "time1", "value")

#Extract max for each vial and add to long datasheet
desmax <-aggregate(des0m[c("value")], by=des0m[c("vial", "treatment")], FUN=max, na.rm=TRUE)
x<-match(des0m$vial, desmax$vial)
des0m$max<-desmax$value[x]

des0m$alive<-(des0m$max-des0m$value)
des0m$time<-as.numeric(as.character(des0m$time1))


dess0<-ddply(des0m, .(treatment,sex,cage), f)
dess0$sample<-0

#Dessication -all other data
des <-read.csv(file= "~/Documents/DrosEED14/Dessication_DrosEEDT1-T4_SMR_170226.csv", stringsAsFactors= TRUE, strip.white = TRUE, na.strings = c("NA","na"), check.names= FALSE)

drop1 <- c("replicate","assayenv", "pheno")
desall = desall[,!(names(desall) %in% drop1)]
View(desall)
##   tubesample cage    treatment sex cagesex  treatmentsex replicate sample
## 1       1one   F1 Non-Seasonal   M     F1M Non-SeasonalM         a    one
## 2       2one   F1 Non-Seasonal   M     F1M Non-SeasonalM         b    one
## 3       3one   F1 Non-Seasonal   M     F1M Non-SeasonalM         c    one
## 4       4one   F1 Non-Seasonal   F     F1F Non-SeasonalF         a    one
## 5       5one   F1 Non-Seasonal   F     F1F Non-SeasonalF         b    one
## 6       6one   F1 Non-Seasonal   F     F1F Non-SeasonalF         c    one
##   cagesexsample vial 4.5 7.5 9.5 11.5 13.5 15.5 17.5 19.5 21.5 23.5 25.5 26.5
## 1        F1Mone F1Ma   0   1   4    9    9   10   10   10   10   10   10   10
## 2        F1Mone F1Mb   0   0   2    8   10   10   10   10   10   10   10   10
## 3        F1Mone F1Mc   0   0   5   10   10   10   10   10   10   10   10   10
## 4        F1Fone F1Fa   0   0   0    2    4    7   11   11   11   11   11   11
## 5        F1Fone F1Fb   0   0   0    2    3    7    7   10   10   10   10   10
## 6        F1Fone F1Fc   0   0   1    2    2    6    9   10   10   10   10   10
##   27.5 28.5 29.5 30.5 31.5 32.5 33.5 34.5 35.5 36.5 37.5 39.5 41.5
## 1   10   10   10   10   10   10   10   10   10   10   10   10   10
## 2   10   10   10   10   10   10   10   10   10   10   10   10   10
## 3   10   10   10   10   10   10   10   10   10   10   10   10   10
## 4   11   11   11   11   11   11   11   11   11   11   11   11   11
## 5   10   10   10   10   10   10   10   10   10   10   10   10   10
## 6   10   10   10   10   10   10   10   10   10   10   10   10   10
des1<-drop.levels(subset(desall, sample=="one"))
View(des1)
des1<-melt(des1)
## Using tubesample, cage, treatment, sex, cagesex, treatmentsex, replicate, sample, cagesexsample, vial as id variables
#Extract max for each vial and add to long datasheet
desmax <-aggregate(des1[c("value")], by=des1[c("vial", "treatment")], FUN=max, na.rm=TRUE)
x<-match(des1$vial, desmax$vial)
des1$max<-desmax$value[x]

des1$alive<-(des1$max-des1$value)
des1$time<-as.numeric(as.character(des1$variable))


dess1<-ddply(des1, .(treatment,sex,cage), f)
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
## Warning: glm.fit: algorithm did not converge
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
dess1$sample<-1

tapply(dess1$V1, list(dess1$treatment, dess1$sex), FUN=mean)
##                     F         M
## Non-Seasonal 14.97600 10.489523
## Seasonal     15.22468  9.639184

#Des2
des2<-drop.levels(subset(desall, sample=="two"))
des2<-melt(des2)
View(des2)
## Using tubesample, cage, treatment, sex, cagesex, treatmentsex, replicate, sample, cagesexsample, vial as id variables
#Extract max for each vial and add to long datasheet
desmax <-aggregate(des2[c("value")], by=des2[c("vial", "treatment")], FUN=max, na.rm=TRUE)
x<-match(des2$vial, desmax$vial)
des2$max<-desmax$value[x]

des2$alive<-(des2$max-des2$value)
des2$time<-as.numeric(as.character(des2$variable))


dess2<-ddply(des2, .(treatment,sex,cage), f)
dess2$sample<-2

tapply(dess2$V1, list(dess2$treatment, dess2$sex), FUN=mean)
##                     F        M
## Non-Seasonal 20.57238 12.25192
## Seasonal     21.96749 12.69166
#Des3
des3<-drop.levels(subset(des, sample=="third"))

des3<-melt(des3)
## Using tubesample, cage, treatment, sex, cagesex, treatmentsex, replicate, sample, cagesexsample, vial as id variables
#Extract max for each vial and add to long datasheet
desmax <-aggregate(des3[c("value")], by=des3[c("vial", "treatment")], FUN=max, na.rm=TRUE)
x<-match(des3$vial, desmax$vial)
des3$max<-desmax$value[x]

des3$alive<-(des3$max-des3$value)
des3$time<-as.numeric(as.character(des3$variable))


dess3<-ddply(des3, .(treatment,sex,cage), f)
dess3$sample<-3

tapply(dess3$V1, list(dess3$treatment, dess3$sex), FUN=mean)
##                     F        M
## Non-Seasonal 20.63213 12.34066
## Seasonal     22.05269 12.76611
dess3F<-drop.levels(subset(dess3, sex=="F"))

z2<-lme(V1 ~ treatment, random = ~1|cage, data=dess3F)

#Des4
des4<-drop.levels(subset(desall, sample=="four"))
des4<-melt(des4)
## Using tubesample, cage, treatment, sex, cagesex, treatmentsex, replicate, sample, cagesexsample, vial as id variables
#Extract max for each vial and add to long datasheet
desmax <-aggregate(des4[c("value")], by=des4[c("vial", "treatment")], FUN=max, na.rm=TRUE)
x<-match(des4$vial, desmax$vial)
des4$max<-desmax$value[x]

des4$alive<-(des4$max-des4$value)
des4$time<-as.numeric(as.character(des4$variable))
View(des2)

dess4<-ddply(des4, .(treatment,sex,cage), f)
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
dess4$sample<-4

tapply(dess4$V1, list(dess4$treatment, dess4$sex), FUN=mean)
##                     F        M
## Non-Seasonal 18.06172 9.715323
## Seasonal     16.96633 9.324134
#Graphs for LD50s for dessication
dess50<-rbind(dess0, dess2, dess4)
View(dess50)
#write.csv(dess50, file="EED14_DessicationLD50_Data_180416.csv")

cdat <- ddply(dess50, c("treatment", "sex", "sample"), summarise,
              N    = sum(!is.na(V1)),
              mean = mean(V1, na.rm=TRUE),
              sd   = sd(V1, na.rm=TRUE),
              se   = sd / sqrt(N)
)


cdat$treatmentsex<-paste(cdat$treatment, cdat$sex)
cdat$treatment<-revalue(cdat$treatment, c("Seasonal"="E", "Non-Seasonal"="N"))

cdatf<-drop.levels(subset(cdat, sex=="F"))
cdat$sample<-as.factor(cdat$sample)

#df_$month<-factor(df_za$month, levels=c("may", "june","july", "august", "october", "february"))

g1 <- ggplot(cdat, aes(x=sample, y=mean, group=treatmentsex))+
  geom_point(aes(shape=sex, color=treatment),size=4)+
  geom_line(aes(color=treatment, linetype=sex))+
  geom_linerange(aes(ymin=mean-se,ymax=mean+se, color=treatment), size=1)+
  ylab("Dessication (hrs) - 2020")+
  xlab("Sample Point")+
  scale_x_discrete(labels=c("Found", "Summer", "Fall"))

g1
## Warning: Removed 2 rows containing missing values (geom_segment).


#Analysis for phenotypic assays that include time

#
#
#
dess0f=subset(dess0,dess0$sex =='F')
dess2f=subset(dess2,dess2$sex == 'F')
dess4f=subset(dess4,dess4$sex == 'F')
View(dess4f)
dess0m=subset(dess0,dess0$sex =='M')
dess2m=subset(dess2,dess2$sex == 'M')
dess4m=subset(dess4,dess4$sex == 'M')
dess50<-rbind(dess0, dess2, dess4)
dess50f<-rbind(dess0f, dess2f, dess4f)
dess50m<-rbind(dess0m, dess2m, dess4m)

dess50$sample<-factor(dess50$sample, levels=c("0", "2", "4"))
dess50f$sample<-factor(dess50f$sample, levels=c("0", "2", "4"))
dess50m$sample<-factor(dess50m$sample, levels=c("0","2","4"))

Dess50<-drop.levels(subset(dess50, cage!='founder' ))
Dess50f<-drop.levels(subset(dess50f, cage!='founder'))
Dess50m<-drop.levels(subset(dess50m, cage!='founder'))

View(Dess50)

x1<-lme(V1 ~ sample*treatment, random=~1|cage/sample, data=Dess50m)
anova(x1)
##                  numDF denDF   F-value p-value
## (Intercept)          1    80 2024.8190  <.0001
## sample               3    54   11.5336  <.0001
## treatment            1    18    0.0894  0.7684
## sample:treatment     3    54    0.4095  0.7468