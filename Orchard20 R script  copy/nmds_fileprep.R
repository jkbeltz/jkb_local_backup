all=NMDS_DF_ALLPHENOS
all$cage.treatment.pheno <- paste(all$Cage,all$Cage.Treatment,all$Pheno.Treatment)
view(all)


## make summary data frame with only APPLE, LB, AT, SAC 024
all.sum=all %>%
  group_by(Label.sum) %>%
  #filter(!(Cage.Treatment == "B"))%>%
  #filter(!any(Timepoint > 1 & Pheno.Treatment == "A")) %>% #Removing apple treatments from bloomington cages
  #filter(!(Timepoint > 1 & Pheno.Treatment == "B")) %>% #Removing bloomington treatments from apple cages
  #filter(!(Pheno.Treatment=="Sterilized"))%>% ##Remove any sterilized phenotyped samples
  #filter(!(Timepoint== 1))%>% ## remove timepoint 1
 #filter(!(Timepoint== 3))%>% ## remove timepoint 3
 #filter(!(Timepoint== 5))%>% ##remove timepoint 5
  summarise(
            Timepoint=as.factor(Timepoint),
            Cage=as.factor(Cage),
            CageTreatment=as.factor(Cage.Treatment),
            PhenoTreatment=as.factor(Pheno.Treatment),
            LD=mean(LD.Mean,na.rm = TRUE),
            LDMale=mean(MLD.Mean,na.rm = TRUE),
            LDFemale=mean(FLD.Mean,na.rm = TRUE),
            Viability=mean(Viability,na.rm = TRUE),
            Starv=mean(Starv.Mean,na.rm = TRUE),
            StarvMale=mean(MStarv.Mean,na.rm = TRUE),
            StarvFemale=mean(FStarv.Mean,na.rm =TRUE),
            Fecundity=mean(F.1DayAVG,na.rm = TRUE),
            DryWeight=mean(`dry weight`,na.rm = TRUE),
            LipidWeight=mean(`lipid weight`,na.rm = TRUE),
            PercentLipid=mean(`%lipid`, na.rm = TRUE),
            Pigmentation=mean(pigmean, na.rm = TRUE)
            
  )
view(all.sum)
all.sum=unique(all.sum)

Fig3DF=subset(all.sum, all.sum$PhenoTreatment=="SAC")
Fig4DF=subset(all.sum, all.sum$PhenoTreatment=="Sterilized")
View(Fig3DF)
View(Fig4DF)

write.csv(Fig3DF, "file:///Users/jackbeltz/Desktop//Fig3DF.csv", row.names=FALSE)
write.csv(Fig4DF, "file:///Users/jackbeltz/Desktop//Fig4DF.csv", row.names=FALSE)

