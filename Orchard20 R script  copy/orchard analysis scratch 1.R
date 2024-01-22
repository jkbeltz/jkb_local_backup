library(ggplot2)
all= Orchard_2020_Phenotype_Master_This
library(tidyverse)
library(dplyr)
all=as_tibble(all)
all.L=subset(all, all$Generation %in% c("F2", "F3"))
all.jack=subset(all.L,all.L$Cage.Treatment %in% c("A","LB","AT"))
all.jack=subset(all.jack,all.jack$Timepoint %in% c(0,2,4))
all.jack=subset(all.jack,all.jack$Pheno.Treatment=="SAC")
View(all.jack)

all.paul=


View(all)
fecund=all %>%
  filter(!any(Assay.ENV== F)) %>% ##remove any field assays
  filter(!any(Timepoint < 1 & Cage.Treatment == "B")) %>%##removing any bloomington cages 
  filter(!any(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  filter(!any(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  filter(!any(Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  filter(!any(Timepoint== 1))%>% ## remove timepoint 1
  filter(!any(Timepoint== 3))%>% ## temove timepoint 3
  filter(!any(Timepoint==5))%>% ##remove timepoint 5
  summarise(treatment=as.factor(Cage.Treatment),
            cage=as.factor(Cage),
            sample=as.factor(Timepoint)
  )
View(fecund)


  