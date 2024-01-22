library(ggplot2)
library(dplyr)
library(tidyr)
##set up DF 
df2 = df_stats_LOA_Meta
View(df2)
#put the matched and target shifts into a single column called 'type'
df = df %>% gather(median.matched, median.target, key = type, value = shift)
#make sure no rows have been duplicated and na's don't exist
df = na.omit(df)
df = unique(df)

### DF Summary table
df.sum=df %>%
  group_by(Test.Seg,Test.Treat,Test.Cage,type) %>%
  dplyr::summarise(N=n(),
            shiftmean=mean(shift,na.rm = TRUE),
            shiftsd=sd(shift,na.rm = TRUE),
            shiftse=shiftsd/sqrt(N)
            
  )
df.sum=unique(df.sum)

###Mark summary viz

df.2_5 = df %>% filter(Test.Seg == '2_5')
df.2_4 = df %>% filter(Test.Seg == '2_4')

ggplot(df.2_5, aes(x = Test.Treat, y = shift, shape = type, colour = type)) +
  geom_jitter(size = 2) +
  scale_shape_manual(values = c(4,19)) +
  scale_color_manual(values = c('light grey', 'orange red')) +
  theme_bw(base_size = 15) +
  #ylim(-0.02, 0.04) +
  theme(legend.position = 'none')

 ####sum of cage
df.sum=df %>%
  group_by(Test.Seg,Test.Treat,Test.Cage,type) %>%
  dplyr::summarise(N=n(),
                   shiftmean=mean(shift,na.rm = TRUE),
                   shiftmedian=median(shift,na.rm = TRUE),
                   shiftsd=sd(shift,na.rm = TRUE),
                   shiftse=shiftsd/sqrt(N)
                   
  )
df.sum=unique(df.sum)

####plot with error bars 
df.sum.2_5 = df.sum %>% filter(Test.Seg == '2_5')
### missing AT 3 and 7(Test cage) 
df.sum.2_4 = df.sum %>% filter(Test.Seg == '2_4')
###missing AT 7 (Test Cage)



ggplot(df.sum.2_5, aes(x = Test.Treat, y = shiftmedian, shape = type, colour = type)) +
  geom_pointrange(aes(ymin=shiftmedian-shiftse,ymax=shiftmedian+shiftse),
                    position=position_jitter(width=.25), 
                    linetype='solid',
                    alpha=.5) +
  scale_shape_manual(values = c(4,19)) +
  scale_color_manual(values = c('light grey', 'orange red')) +
  theme_bw(base_size = 15) +
  #ylim(-0.02, 0.04) +
  theme(legend.position = 'none')

####Sum of treatment 
df.sumtreat=df %>%
  group_by(Test.Seg,Test.Treat,type) %>%
  dplyr::summarise(N=n(),
                   shiftmean=mean(shift,na.rm = TRUE),
                   shiftmedian=median(shift,na.rm = TRUE),
                   shiftsd=sd(shift,na.rm = TRUE),
                   shiftse=shiftsd/sqrt(N)
                   
  )
df.sumtreat=unique(df.sumtreat)
#View(df.sumtreat)
####plot with error bars 
df.sumtreat.2_5 = df.sumtreat %>% filter(Test.Seg == '2_5')
### missing AT 3 and 7(Test cage) 
df.sumtreat.2_4 = df.sumtreat %>% filter(Test.Seg == '2_4')
###missing AT 7 (Test Cage)



ggplot(df.sumtreat.2_5, aes(x = Test.Treat, y = shiftmedian, shape = type, colour = type)) +
  geom_pointrange(aes(ymin=shiftmedian-shiftse,ymax=shiftmedian+shiftse),
                  linetype='solid',
                  size=2) +
  geom_pointrange(aes(ymin=shiftmedian-shiftse,ymax=shiftmedian+shiftse),
                  position=position_jitter(width=.2), 
                  linetype='solid',
                  alpha=.3,
                  data = df.sum.2_5) +
  scale_shape_manual(values = c(4,19)) +
  scale_color_manual(values = c('light grey', 'orange red')) +
  theme_bw(base_size = 15) +
  #ylim(-0.02, 0.04) +
  theme(legend.position = 'none')

### sum by LOC
df.sumLOC=df %>%
  group_by(Test.Seg,Test.Treat,LOC,type) %>%
  dplyr::summarise(N=n(),
                   shiftmean=mean(shift,na.rm = TRUE),
                   shiftmedian=median(shift,na.rm = TRUE),
                   shiftsd=sd(shift,na.rm = TRUE),
                   shiftse=shiftsd/sqrt(N)
                   
  )
df.sumLOC=unique(df.sumLOC)
View(df.sumLOC)
####plot with error bars 
df.sumLOC.2_5 = df.sumLOC %>% filter(Test.Seg == '2_5')
### missing AT 3 and 7(Test cage) 
df.sumLOC.2_4 = df.sumLOC %>% filter(Test.Seg == '2_4')
###missing AT 7 (Test Cage)



ggplot(df.sumLOC.2_5, aes(x = Test.Treat, y = shiftmedian, shape = type, colour = type)) +
  geom_pointrange(aes(ymin=shiftmedian-shiftse,ymax=shiftmedian+shiftse),
                  position=position_jitter(width=.25), 
                  linetype='solid') +
  scale_shape_manual(values = c(4,19)) +
  scale_color_manual(values = c('light grey', 'orange red')) +
  theme_bw(base_size = 15) +
  #ylim(-0.02, 0.04) +
  theme(legend.position = 'none')



















####NEW DF####


df=df_shifts_LOA_Meta_FDR1_Shift01
df$Type.Treat=paste(df$Type, df$Test.Treatment)

####SUBSETS
df.mAT=as.numeric(unlist(subset(df[,11], df$Type.Treat =="matched AT")))
df.mLB=as.numeric(unlist(subset(df[,11], df$Type.Treat =="matched LB")))
df.mMA=as.numeric(unlist(subset(df[,11], df$Type.Treat == c("matched AT","matched LB"))))

#df <- subset(df , !df$Type == "matched" | !df$Test.Treatment %in% c("AT", "LB")) ##remove AT and LB matched
df.mA=as.numeric(unlist(subset(df[,11], df$Type.Treat =="matched A")))
df.tA=as.numeric(unlist(subset(df[,11], df$Type.Treat =="target A")))
df.tAT=as.numeric(unlist(subset(df[,11], df$Type.Treat =="target AT")))
df.tLB=as.numeric(unlist(subset(df[,11], df$Type.Treat =="target LB")))
df.tMA=as.numeric(unlist(subset(df[,11], df$Type.Treat == c("target AT","target LB"))))

df2L <- subset (df , df$chrom == "2L")
df2L.mA=as.numeric(unlist(subset(df2L[,11], df2L$Type.Treat =="matched A")))
df2L.tA=as.numeric(unlist(subset(df2L[,11], df2L$Type.Treat =="target A")))
df2L.tAT=as.numeric(unlist(subset(df2L[,11], df2L$Type.Treat =="target AT")))
df2L.tLB=as.numeric(unlist(subset(df2L[,11], df2L$Type.Treat =="target LB")))

df2R <- subset (df , df$chrom == "2R")
df2R.mA=as.numeric(unlist(subset(df2R[,11], df2R$Type.Treat =="matched A")))
df2R.tA=as.numeric(unlist(subset(df2R[,11], df2R$Type.Treat =="target A")))
df2R.tAT=as.numeric(unlist(subset(df2R[,11], df2R$Type.Treat =="target AT")))
df2R.tLB=as.numeric(unlist(subset(df2R[,11], df2R$Type.Treat =="target LB")))

df3L <- subset (df , df$chrom == "3L")
df3L.mA=as.numeric(unlist(subset(df3L[,11], df3L$Type.Treat =="matched A")))
df3L.tA=as.numeric(unlist(subset(df3L[,11], df3L$Type.Treat =="target A")))
df3L.tAT=as.numeric(unlist(subset(df3L[,11], df3L$Type.Treat =="target AT")))
df3L.tLB=as.numeric(unlist(subset(df3L[,11], df3L$Type.Treat =="target LB")))

df3R <- subset (df , df$chrom == "3R")
df3R.mA=as.numeric(unlist(subset(df3R[,11], df3R$Type.Treat =="matched A")))
df3R.tA=as.numeric(unlist(subset(df3R[,11], df3R$Type.Treat =="target A")))
df3R.tAT=as.numeric(unlist(subset(df3R[,11], df3R$Type.Treat =="target AT")))
df3R.tLB=as.numeric(unlist(subset(df3R[,11], df3R$Type.Treat =="target LB")))

dfX <- subset (df , df$chrom == "X")
dfX.mA=as.numeric(unlist(subset(dfX[,11], dfX$Type.Treat =="matched A")))
dfX.tA=as.numeric(unlist(subset(dfX[,11], dfX$Type.Treat =="target A")))
dfX.tAT=as.numeric(unlist(subset(dfX[,11], dfX$Type.Treat =="target AT")))
dfX.tLB=as.numeric(unlist(subset(dfX[,11], dfX$Type.Treat =="target LB")))


####KS TEST####

##targets##
a=ks.test(df.mA, df.mAT)
a
##df.mAT vs. df.mLB - D = 0.010851, p-value < 2.2e-16
##df.mA vs. df.mLB - D = 0.019227, p-value < 2.2e-16
##df.mA and df.mAT - D = 0.029612, p-value < 2.2e-16

## A vs MA (AT+LB)
a=ks.test(df.tA, df.tMA)
a






##all chromosomes - everything very sig
a=ks.test(df.tAT, df.tLB)
length(df.LB)
##all sig chromo combined
## df.mA vs. df.tA  --- D = 0.16588, p-value < 2.2e-16
## df.tA and df.tAT --- D = 0.03582 - 0.029612 = 0.006208, p-value < 2.2e-16
## df.tA and df.tLB --- D = 0.048696 - 0.019227 = 0.029469, p-value < 2.2e-16
## df.tAT and df.tLB --- D = 0.041505 - 0.010851 = 0.030654 , p-value < 2.2e-16

## 2L 
a=ks.test(df2L.tAT, df2L.tLB)

##mA vs tA  - D = 0.22522, p-value < 2.2e-16
##tA vs tAT - D = 0.028348, p-value < 2.2e-16
##tA vs tLB - D = 0.079198, p-value < 2.2e-16
##tAT vs tLB - D = 0.054675, p-value < 2.2e-16

## 2R
a=ks.test(df2R.tAT, df2R.tLB)

##mA vs tA  - D = 0.078534, p-value < 2.2e-16
##tA vs tAT - D = 0.11508, p-value < 2.2e-16
##tA vs tLB - D = 0.089669, p-value < 2.2e-16
##tAT vs tLB - D = 0.028798, p-value < 2.2e-16

##3L
a=ks.test(df3L.tAT, df3L.tLB)

##mA vs tA  - D = 0.10632, p-value < 2.2e-16
##tA vs tAT - D = 0.073108, p-value < 2.2e-16
##tA vs tLB - D = 0.11614, p-value < 2.2e-16
##tAT vs tLB - D = 0.060998, p-value < 2.2e-16

##3R
a=ks.test(df3R.tAT, df3R.tLB)

##mA vs tA  - D = 0.2239, p-value < 2.2e-16
##tA vs tAT - D = 0.055012, p-value < 2.2e-16
##tA vs tLB - D = 0.12842, p-value < 2.2e-16
##tAT vs tLB - D = 0.09761, p-value < 2.2e-16

##X
a=ks.test(dfX.tAT, dfX.tLB)
a
##mA vs tA  - D = 0.078926, p-value < 2.2e-16
##tA vs tAT -D = 0.08007, p-value < 2.2e-16
##tA vs tLB - D = 0.077795, p-value < 2.2e-16
##tAT vs tLB - D = 0.036886, p-value < 2.2e-16

df.AATLB=df %>%
  filter(Type.Treat == "matched A"|
         Type.Treat == "target A"|
         Type.Treat == "target AT"|
         Type.Treat == "target LB")


View(df.AATLB)

#### Summary Files 

dfch.sum <- df.AATLB %>%
  group_by(Type.Treat, chrom) %>%
  summarize(mean = mean(TestCageShift),
            median = median(TestCageShift))


df.sum <- df %>%
  group_by(Type.Treat) %>%
  summarize(mean = mean(TestCageShift),
            median = median(TestCageShift))

## By chrom 
df.sum2L= subset(dfch.sum, dfch.sum$chrom =="2L")
df.sum2R= subset(dfch.sum, dfch.sum$chrom =="2R")
df.sum3L= subset(dfch.sum, dfch.sum$chrom =="3L")
df.sum3R= subset(dfch.sum, dfch.sum$chrom =="3R")


##make MA 
library(dplyr)
MAdf=df %>%
  mutate(Type.Treat.MA = case_when(Type.Treat == "matched A" ~ 'matched A',
                                   Type.Treat == "target A" ~ 'target A',
                                   #Type.Treat == "matched AT" ~ 'matched MA',
                                   #Type.Treat == "matched LB" ~ 'matched MA',
                                   Type.Treat == "target LB" ~ 'target MA',
                                   Type.Treat == "target AT" ~ 'target MA',))

MAdf=na.omit(MAdf)


MAdfch.sum=dfch.sum %>%
  mutate(Type.Treat.MA = case_when(Type.Treat == "matched A" ~ 'matched A',
                                   Type.Treat == "target A" ~ 'target A',
                                   #Type.Treat == "matched AT" ~ 'matched MA',
                                   #Type.Treat == "matched LB" ~ 'matched MA',
                                   Type.Treat == "target LB" ~ 'target MA',
                                   Type.Treat == "target AT" ~ 'target MA ',))

View (MAdfch.sum)

##this plot
 r=ggplot(df.AATLB, aes(x=TestCageShift, color = Type.Treat)) +
  geom_density(size=.5) +
  facet_wrap(. ~ chrom) +
  scale_color_manual(values = c("matched A"="black","target A"="#848FA2","target AT"="#CC2D35","target LB"="#058ED9"))+
  geom_vline(data = dfch.sum, aes(xintercept = median, color = Type.Treat ),size=.1)+
  theme_bw()

r + coord_cartesian(xlim=c(-0.2,0.2), ylim = c(0,25)) + theme_classic()
r + coord_cartesian(xlim=c(-0.001,0.04), ylim = c(10,25)) + theme_classic()

##A vs, MA

d=ggplot(MAdf, aes(x=TestCageShift, color = Type.Treat.MA)) +
  geom_density(size=1) +
  facet_wrap(. ~ chrom) +
  scale_color_manual(values = c("matched A"="black","target A"="#848FA2","target MA"="purple"))+
  #geom_vline(data = MAdfch.sum, aes(xintercept = median, color = Type.Treat.MA ),size=.2)+
  theme_classic()

d+ coord_cartesian(xlim=c(-0.2,0.2), ylim = c(0,25)) + theme_classic()
d + coord_cartesian(xlim=c(-0.001,0.04), ylim = c(10,25)) + theme_classic()


#3filled plot
ggplot(df, aes(x=TestCageShift, fill = Type.Treat)) +
  geom_density(size=.5, alpha = .6) +
  facet_wrap(. ~ chrom) +
  scale_fill_manual(values = c("matched A"="black","target A"="#848FA2","target AT"="#CC2D35","target LB"="#058ED9"))+
  geom_vline(data = df.sum, aes(xintercept = median, color = Type.Treat ),size=.1)+
  theme_bw()
