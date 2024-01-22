#### orchard 2020 16s analysis #########################
## started 9/14/21 by jk beltz
      ## post qiime2 analysis ########
install.packages("tidyverse")
library("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

BiocManager::install("dada2", version = "3.17")

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.17")

#!/usr/bin/env Rscript
library(qiime2R)
library(tidyverse)
library(Biostrings)
library(dplyr)
install.packages("dada2")
library(dada2)
library(qiime2R)
library(ggplot2)

### online tutorial https://www.google.com/search?q=qiimer+vs.+qiime2r&oq=qiimer+vs.+qiime2r&aqs=chrome..69i57j0i22i30.11627j0j7&sourceid=chrome&ie=UTF-8
##Reading QZA files#######
##SET WD to HOME/DOCUMENTS/16S_QIIME2
?read_qza
SVs<-read_qza("table.qza")
View(SVs)
OTUTABLE=SVs$data


View(OTUTABLE)
SVs$data[1:5,1:5]
SVs$uuid
SVs$type
SVs$contents
print_provenance(SVs)

denoise_stats<-read_qza("Orchard20/QIIME_output/denoise_fwd_0-240_rev_0-240/denoise_stats.qza")
rep_seqs=read_qza("Orchard20/QIIME_output/denoise_fwd_0-240_rev_0-240/representative-seqs.qza")
unrooted_artifact=read_qza("Orchard20/QIIME_output/denoise_fwd_0-240_rev_0-240/unrooted-tree.qza")
rooted_artifact=read_qza("Orchard20/QIIME_output/denoise_fwd_0-240_rev_0-240/rooted-tree.qza")
taxonomy=read_qza("correcttaxonomy.qza")

head(taxonomy$data)
taxon<-parse_taxonomy(taxonomy$data) #### UNIQUE COLUMN FOR EACH TAXON LEVEL
View(taxon)

## READ .BIOM #########
?read_q2biom
taxontable=read_q2biom("table-taxonomy2.biom") ##as.biom
View(taxontable)
TaxonTable=read.csv(file='orchard20otutable.csv', header = TRUE) ##as.csv
colnames(TaxonTable)[1] <- "Feature.ID"### changing names of columns for parse function
colnames(TaxonTable)[244] <- "Taxon"
TaxonTable$Taxon
TaxonTable<-parse_taxonomy(TaxonTable$taxonomy) ## parse still didnt work so merged raw .qzas based on row name

otu.taxon<- merge(OTUTABLE, taxon, by=0, all=TRUE) ## successful merge, names double checked
ncol(OTUTABLE_nw)
length(metadata)

View(otu.taxon)
metadata=read.csv(file = "orchard2016smetadata.csv", header = TRUE)
colnames(metadata)[1]= "SampleID"
View(metadata)
## loading in qiime2 metric files ###########

evenvess=scan("metrics/evenness_vector.qza") ##error 
shannon=read_qza("metrics/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID")
observed_features=read_qza("metrics/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID")
faith_pd=read_qza("metrics/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID")
rarefied=read_qza("metrics/rarefied_table.qza")
## cant rename column b/c sample names are variables 
uu_distance_matrix=read_qza("metrics/unweighted_unifrac_distance_matrix.qza")
#View(uu_distance_matrix$data)
bray_distance_matrix=read_qza("metrics/bray_curtis_distance_matrix.qza")
bray_curtis_pcoa=read_qza("metrics/bray_curtis_pcoa_results.qza")
jaccard_matrix=read_qza("metrics/jaccard_distance_matrix.qza")
jaccard_pcoa=read_qza("metrics/jaccard_pcoa_results.qza")
uu_pcoa=read_qza("metrics/unweighted_unifrac_pcoa_results.qza")
wu_distance_matrix=read_qza("metrics/weighted_unifrac_distance_matrix.qza")
wu_pcoa=read_qza("metrics/weighted_unifrac_pcoa_results.qza")


########ALPHA DIVERSITY#############

View(shannon$data)
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))
### 14 unaccounted files in the shannon file

write.csv(shannon, file = "Orchard20RWD/shannon.csv") ##compared them directly in numbers, created file to do inventory of whihc resultys are missing from each metric  
write.csv(metadata, file ="Orchard20RWD/metadata.csv")


observed_features<-observed_features$data %>% rownames_to_column("SampleID")
gplots::venn(list(metadata=metadata$SampleID, observed_features=observed_features$SampleID))

write.csv(observed_features, file ="Orchard20RWD/observed_features.csv")
#same missing 14 values - confirmed these value have no detected reads, meaningful result. 

faith_pd<-faith_pd$data %>% rownames_to_column("SampleID")

######removed wolbachia in qiime, see "useful Qiime Commands"############
##only removed 1 sample from analysis and it was a control!

observed_features_nw=read_qza("metrics-no-wolbachia/observed_features_vector.qza")
observed_features_nw<-observed_features_nw$data %>% rownames_to_column("SampleID")
write.csv(observed_features_nw, "Orchard20RWD/observed_features_nw.csv")

##reload all data with wolbachia samples removed

evenvess_nw=scan("metrics-no-wolbachia/evenness_vector.qza") ##error 
shannon_nw=read_qza("metrics-no-wolbachia/shannon_vector.qza")
shannon_nw<-shannon_nw$data %>% rownames_to_column("SampleID")
#View(shannon_nw)
observed_features_nw=read_qza("metrics-no-wolbachia/observed_features_vector.qza")
observed_features_nw<-observed_features_nw$data %>% rownames_to_column("SampleID")
faith_pd_nw=read_qza("metrics-no-wolbachia/faith_pd_vector.qza")
faith_pd_nw<-faith_pd_nw$data %>% rownames_to_column("SampleID")
rarefied_nw=read_qza("metrics-no-wolbachia/rarefied_table.qza")

## cant rename column b/c sample names are variables 
uu_distance_matrix_nw=read_qza("metrics-no-wolbachia/unweighted_unifrac_distance_matrix.qza")
uu_distance_matrix_nw=uu_distance_matrix_nw$data
bray_distance_matrix_nw=read_qza("metrics-no-wolbachia/bray_curtis_distance_matrix.qza")
bray_distance_matrix_nw=bray_distance_matrix_nw$data
bray_curtis_pcoa_nw=read_qza("metrics-no-wolbachia/bray_curtis_pcoa_results.qza")
bray_curtis_pcoa_nw=bray_curtis_pcoa_nw$data
jaccard_matrix_nw=read_qza("metrics-no-wolbachia/jaccard_distance_matrix.qza")
jaccard_matrix_nw=jaccard_matrix_nw$data
jaccard_pcoa_nw=read_qza("metrics-no-wolbachia/jaccard_pcoa_results.qza")
jaccard_pcoa_nw=jaccard_pcoa_nw$data
uu_pcoa_nw=read_qza("metrics-no-wolbachia/unweighted_unifrac_pcoa_results.qza")
uu_pcoa_nw=uu_pcoa_nw$data
wu_distance_matrix_nw=read_qza("metrics-no-wolbachia/weighted_unifrac_distance_matrix.qza")
wu_distance_matrix_nw=wu_distance_matrix_nw$data
wu_pcoa_nw=read_qza("metrics-no-wolbachia/weighted_unifrac_pcoa_results.qza")
#View(wu_pcoa_nw)

####reread otutable 

SVs_nw=read_qza("table-no-wolbachia-exact.qza")
OTUTABLE_nw=SVs_nw$data
typeof(OTUTABLE_nw)
#View(taxon)

#remerge old tacxon table with new otu values
otu.taxon_nw<- merge(taxon, OTUTABLE_nw, by ='row.names', all=TRUE) ## successful merge, names double checked

#View(otu.taxon_nw)
######metadata work#############
metadata=read.csv(file = "orchard2016smetadata.csv", header = TRUE)
colnames(metadata)[1]= "SampleID"
#View(metadata)
metadata[,20]=tolower(metadata[,20])
metadata[,22]=tolower(metadata[,22])

### add alpha diveristy metrics to meta###################
a=metadata %>% 
  merge(observed_features_nw, all.y=TRUE) 
b=a %>% 
  merge(faith_pd_nw, all.y=TRUE) 
allmeta=b %>% 
  merge(shannon_nw, all.y=TRUE) 
View (allmeta)
alphadata = allmeta
#alphadata <- allmeta[ -c(2:7, 9:13, 23:36, 39:49) ]
#alphadata[101,5]="Axenic_LB"
View(alphadata)

alphadatano5<-alphadata[(alphadata$time_point == "T2") | (alphadata$time_point == "T4"),]

###subsetting####
APBloomalpha<-alphadatano5[(alphadatano5$cage_treatment == "Apple") | (alphadatano5$cage_treatment == "Bloom"),]
APATLBalpha<-alphadatano5[(alphadatano5$cage_treatment == "Apple") | (alphadatano5$cage_treatment == "Apple_AT") | (alphadatano5$cage_treatment == "Apple_LB"),]

#cagefood
APBloomfoodalpha=APBloomalpha[(APBloomalpha$sample_type_detail== "cage food"),]
APATLBfoodalpha=APATLBalpha[(APATLBalpha$sample_type_detail== "cage food"),]

#wild flies 
APBloomwildalpha=APBloomalpha[(APBloomalpha$sample_type_detail== "flies"),]
APATLBwildalpha=APATLBalpha[(APATLBalpha$sample_type_detail== "flies"),]

#F1s
APBloomF1salpha=APBloomalpha[(APBloomalpha$sample_type_detail== "dry stored f1s"),]
APATLBF1salpha=APATLBalpha[(APATLBalpha$sample_type_detail== "dry stored f1s"),]


alphacageenv= subset(alphadata, alphadata$SampleType=="Environmental control")
alphaallflies=subset(alphadata, alphadata$SampleType=="Whole fly")
alphawildflies=subset(alphadata, alphadata$sample_type_detail== "flies" )
#alphawildF1s=subset(alphadata, alphadata$sample_type_detail= c("flies", "dry stored f1s"))
alphawildfood<-alphadata[(alphadata$sample_type_detail== "flies") | (alphadata$sample_type_detail == "cage food"),]
alphawildF1s <-alphadata[(alphadata$sample_type_detail == "flies") | (alphadata$sample_type_detail == "dry stored f1s"),]
alphawildfoodf1s <-alphadata[(alphadata$sample_type_detail == "flies") | (alphadata$sample_type_detail =="dry stored f1s") | (alphadata$sample_type_detail== "cage food"),]
alphaf2s=subset(alphadata, alphadata$sample_type_detail == "dry stored f2s")

alphaallflies$cage_treatment= as.factor(alphaallflies$cage_treatment)
alphaallflies$time_point=as.factor(alphaallflies$time_point)
alphaallflies$sample_type_detail=as.factor(alphaallflies$sample_type_detail)
###line plotting#####

ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


ggplot(data=alphacageenv, aes(x=as.factor('time_point'), y='observed_features', group_by='cage_treatment')) +
  geom_line()+
  geom_point()

p <- ggplot(alphacageenv, aes(x =time_point, y = shannon_entropy, color = cage_treatment))
# Change line types and point shapes by groups
p + geom_line(aes(linetype = cage_treatment)) +
  geom_point(aes(shape = sample_type_detail))

alphawildflies$cage_treatment= as.factor(alphawildflies$cage_treatment)
alphawildflies$time_point=as.factor(alphawildflies$time_point)
alphawildflies$sample_type_detail=as.factor(alphawildflies$sample_type_detail)
alphawildflies$study_group=as.factor(alphawildflies$study_group)

#View(APATLBF1salpha)

#### best plot, just change dataframe #################
APATLBwildalpha%>%
  filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=shannon_entropy, group=interaction(as.factor(cage_treatment), as.factor(collection_type)), color=as.factor(cage_treatment), shape=as.factor(sample_type_detail))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=4, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,5)+
  xlab("Time Point") +
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Apple_LB", "Apple_AT"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations"),
                      values = c("#848FA2","#058ED9","#CC2D35"))+
  theme_classic()  # try other themes like theme_bw() or theme_classic()
  #scale_color_viridis_d(name="Sample Type") # use different color scale which is color blind friendly

APBloomF1salpha%>%
  filter (!sample_type_detail == 'cage env') %>%
  ggplot(aes(x=as.factor(time_point), y=shannon_entropy, group=interaction(as.factor(cage_treatment), as.factor(collection_type)), color=as.factor(cage_treatment), shape=as.factor(sample_type_detail))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  #geom_jitter(shape=0, width=0.2, height=0)+
  stat_summary(geom="point", fun.data=mean_se, size=4, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,5)+
  xlab("Time Point") +
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple", "Bloom"),
                      labels = c("Apple","Bloomington"),
                      values = c("#CC2D35","#000000"))+
  theme_classic()  # try other themes like theme_bw() or theme_classic()
#scale_color_viridis_d(name="Sample Type") # use different color scale which is color blind friendly






#ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

###BAR wit jitter
alphawildfoodf1s %>%
  filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(x=time_point, y=faith_pd , fill= sample_type_detail)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black", position = position_dodge()) + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  facet_grid(~ cage_treatment) + # create a panel for each body site
  xlab("Sample Timepoints") +
  ylab("Observed Features Count") +
  theme_q2r() 
  



ggsave("../../../images/Shannon_by_abx.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches




#########################PCAs#############################################################################
metadata<-read_q2metadata("metric/sample-metadata.tsv")
uwunifrac<-read_qza("metrics/unweighted_unifrac_pcoa_results.qza")
shannon<-read_qza("metrics/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

uualphameta=uu_pcoa_nw$Vectors %>%
  select(SampleID, PC1, PC2, PC3, PC4) %>%
  left_join(alphadata)

wualphameta=wu_pcoa_nw$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3, PC4) %>%
  left_join(alphadata)

View(uualphameta)
View(alphadata)
##subsets #########
##Wild flies + cage food 
uualphawildcage=uualphameta[(uualphameta$sample_type_detail == "flies") | (uualphameta$sample_type_detail == "cage food"),]

uualphawildcage2=uualphawildcage[(uualphawildcage$time_point == "T2"),] 
uualphawildcage4=uualphawildcage[(uualphawildcage$time_point == "T4"),] 


uualphawildcageAB=uualphawildcage[(uualphawildcage$cage_treatment == "Apple") | (uualphawildcage$cage_treatment == "Bloom"),]
uualphawildcageATLB=uualphawildcage[(uualphawildcage$cage_treatment == "Apple_AT") | (uualphawildcage$cage_treatment == "Apple_LB"),]
uualphawildcageA=uualphawildcage[(uualphawildcage$cage_treatment == "Apple"),]
uualphawildcageB=uualphawildcage[(uualphawildcage$cage_treatment == "Bloom"),]
uualphawildcageAT=uualphawildcage[(uualphawildcage$cage_treatment == "Apple_AT"),]
uualphawildcageLB=uualphawildcage[(uualphawildcage$cage_treatment == "Apple_LB"),]

wualphawildcage=wualphameta[(wualphameta$sample_type_detail == "flies") | (wualphameta$sample_type_detail == "cage food"),]

wualphawildcage2=wualphawildcage[(wualphawildcage$time_point == "T2"),] 
wualphawildcage4=wualphawildcage[(wualphawildcage$time_point == "T4"),]

wualphawildcageAB=wualphawildcage[(wualphawildcage$cage_treatment == "Apple") | (wualphawildcage$cage_treatment == "Bloom"),]
wualphawildcageATLB=wualphawildcage[(wualphawildcage$cage_treatment == "Apple_AT") | (wualphawildcage$cage_treatment == "Apple_LB"),]
wualphawildcageA=wualphawildcage[(wualphawildcage$cage_treatment == "Apple"),]
wualphawildcageB=wualphawildcage[(wualphawildcage$cage_treatment == "Bloom"),]
wualphawildcageAT=wualphawildcage[(wualphawildcage$cage_treatment == "Apple_AT"),]
wualphawildcageLB=wualphawildcage[(wualphawildcage$cage_treatment == "Apple_LB"),]

## wild flies + f1s
uualphameta[10,9]="base Bloom"
uualphameta[11,9]="base Bloom"
#View(uualphameta)

uualphawildF1s=uualphameta[(uualphameta$sample_type_detail == "flies") | (uualphameta$sample_type_detail == "dry stored f1s"),]

uualphawildF1s2=uualphawildF1s[(uualphawildF1s$time_point == "T2"),]
uualphawildF1s4=uualphawildF1s[(uualphawildF1s$time_point == "T4"),]

uualphawildF1sAB=uualphawildF1s[(uualphawildF1s$cage_treatment == "Apple") | (uualphawildF1s$cage_treatment == "Bloom"),]
uualphawildF1sATLB=uualphawildF1s[(uualphawildF1s$cage_treatment == "Apple_AT") | (uualphawildF1s$cage_treatment == "Apple_LB"),]
uualphawildF1sA=uualphawildF1s[(uualphawildF1s$cage_treatment == "Apple"),]
uualphawildF1sB=uualphawildF1s[(uualphawildF1s$cage_treatment == "Bloom"),]
uualphawildF1sAT=uualphawildF1s[(uualphawildF1s$cage_treatment == "Apple_AT"),]
uualphawildF1sLB=uualphawildF1s[(uualphawildF1s$cage_treatment == "Apple_LB"),]

wualphameta[10,9]="base Bloom"
wualphameta[11,9]="base Bloom"

wualphawildF1s=wualphameta[(wualphameta$sample_type_detail == "flies") | (wualphameta$sample_type_detail == "dry stored f1s"),]

wualphawildF1s2=wualphawildF1s[(wualphawildF1s$time_point == "T2"),]
wualphawildF1s4=wualphawildF1s[(wualphawildF1s$time_point == "T4"),]

wualphawildF1sAB=wualphawildF1s[(wualphawildF1s$cage_treatment == "Apple") | (wualphawildF1s$cage_treatment == "Bloom"),]
wualphawildF1sATLB=wualphawildF1s[(wualphawildF1s$cage_treatment == "Apple_AT") | (wualphawildF1s$cage_treatment == "Apple_LB"),]
wualphawildF1sA=wualphawildF1s[(wualphawildF1s$cage_treatment == "Apple"),]
wualphawildF1sB=wualphawildF1s[(wualphawildF1s$cage_treatment == "Bloom"),]
wualphawildF1sAT=wualphawildF1s[(wualphawildF1s$cage_treatment == "Apple_AT"),]
wualphawildF1sLB=wualphawildF1s[(wualphawildF1s$cage_treatment == "Apple_LB"),]

## wild flies alone
wualphawild=wualphameta[(wualphameta$sample_type_detail == "flies") ,]

wualphawild2=wualphawild[(wualphawild$time_point == "T2"),]
wualphawild4=wualphawild[(wualphawild$time_point == "T4"),]

wualphawildAB=wualphawild[(wualphawild$cage_treatment == "Apple") | (wualphawild$cage_treatment == "Bloom"),]
wualphawildATLB=wualphawild[(wualphawild$cage_treatment == "Apple_AT") | (wualphawild$cage_treatment == "Apple_LB"),]
wualphawildA=wualphawild[(wualphawild$cage_treatment == "Apple"),]
wualphawildB=wualphawild[(wualphawild$cage_treatment == "Bloom"),]
wualphawildAT=wualphawild[(wualphawild$cage_treatment == "Apple_AT"),]
wualphawildLB=wualphawild[(wualphawild$cage_treatment == "Apple_LB"),]

uualphawild=uualphameta[(uualphameta$sample_type_detail == "flies") ,]

uualphawild2=uualphawild[(uualphawild$time_point == "T2"),]
uualphawild4=uualphawild[(uualphawild$time_point == "T4"),]

uualphawildAB=uualphawild[(uualphawild$cage_treatment == "Apple") | (uualphawild$cage_treatment == "Bloom"),]
uualphawildATLB=uualphawild[(uualphawild$cage_treatment == "Apple_AT") | (uualphawild$cage_treatment == "Apple_LB"),]
uualphawildA=uualphawild[(uualphawild$cage_treatment == "Apple"),]
uualphawildB=uualphawild[(uualphawild$cage_treatment == "Bloom"),]
uualphawildAT=uualphawild[(uualphawild$cage_treatment == "Apple_AT"),]
uualphawildLB=uualphawild[(uualphawild$cage_treatment == "Apple_LB"),]

## F1s alone
wualphaF1s=wualphameta[(wualphameta$sample_type_detail == "dry stored f1s") ,]

wualphaF1s2=wualphaF1s[(wualphaF1s$time_point == "T2"),]
wualphaF1s4=wualphaF1s[(wualphaF1s$time_point == "T4"),]

wualphaF1sAB=wualphaF1s[(wualphaF1s$cage_treatment == "Apple") | (wualphaF1s$cage_treatment == "Bloom"),]
wualphaF1sATLB=wualphaF1s[(wualphaF1s$cage_treatment == "Apple_AT") | (wualphaF1s$cage_treatment == "Apple_LB"),]
wualphaF1sA=wualphaF1s[(wualphaF1s$cage_treatment == "Apple"),]
wualphaF1sB=wualphaF1s[(wualphaF1s$cage_treatment == "Bloom"),]
wualphaF1sAT=wualphaF1s[(wualphaF1s$cage_treatment == "Apple_AT"),]
wualphaF1sLB=wualphaF1s[(wualphaF1s$cage_treatment == "Apple_LB"),]

uualphaF1s=uualphameta[(uualphameta$sample_type_detail == "dry stored f1s") ,]

uualphaF1s2=uualphaF1s[(uualphaF1s$time_point == "T2"),]
uualphaF1s4=uualphaF1s[(uualphaF1s$time_point == "T4"),]

uualphaF1sAB=uualphaF1s[(uualphaF1s$cage_treatment == "Apple") | (uualphaF1s$cage_treatment == "Bloom"),]
uualphaF1sATLB=uualphaF1s[(uualphaF1s$cage_treatment == "Apple_AT") | (uualphaF1s$cage_treatment == "Apple_LB"),]
uualphaF1sA=uualphaF1s[(uualphaF1s$cage_treatment == "Apple"),]
uualphaF1sB=uualphaF1s[(uualphaF1s$cage_treatment == "Bloom"),]
uualphaF1sAT=uualphaF1s[(uualphaF1s$cage_treatment == "Apple_AT"),]
uualphaF1sLB=uualphaF1s[(uualphaF1s$cage_treatment == "Apple_LB"),]

## F2s alone
wualphaF2s=wualphameta[(wualphameta$sample_type_detail == "dry stored f2s") ,]

wualphaF2s2=wualphaF2s[(wualphaF2s$time_point == "T2"),]
wualphaF2s4=wualphaF2s[(wualphaF2s$time_point == "T4"),]

wualphaF2sATLB=wualphaF2s[(wualphaF2s$cage_treatment == "Axenic_AT") | (wualphaF2s$cage_treatment == "Axenic_LB"),]
wualphaF2sAT=wualphaF2s[(wualphaF2s$cage_treatment == "Axenic_AT"),]
wualphaF2sLB=wualphaF2s[(wualphaF2s$cage_treatment == "Axenic_LB"),]

uualphaF2s=uualphameta[(uualphameta$sample_type_detail == "dry stored f2s") ,]

uualphaF2s2=uualphaF2s[(uualphaF2s$time_point == "T2"),]
uualphaF2s4=uualphaF2s[(uualphaF2s$time_point == "T4"),]

uualphaF2sATLB=uualphaF2s[(uualphaF2s$cage_treatment == "Axenic_AT") | (uualphaF2s$cage_treatment == "Axenic_LB"),]

uualphaF2sAT=uualphaF2s[(uualphaF2s$cage_treatment == "Axenic_AT"),]
uualphaF2sLB=uualphaF2s[(uualphaF2s$cage_treatment == "Axenic_LB"),]

## cage food + cage env
wualphaenvfood=wualphameta[(wualphameta$sample_type_detail == "cage food") | (wualphameta$sample_type_detail == "cage env"),]

wualphaenvfood2=wualphaenvfood[(wualphaenvfood$time_point == "T2"),]
wualphaenvfood4=wualphaenvfood[(wualphaenvfood$time_point == "T4"),]

wualphaenvfoodAB=wualphaenvfood[(wualphaenvfood$cage_treatment == "Apple") | (wualphaenvfood$cage_treatment == "Bloom"),]
wualphaenvfoodATLB=wualphaenvfood[(wualphaenvfood$cage_treatment == "Apple_AT") | (wualphaenvfood$cage_treatment == "Apple_LB"),]
wualphaenvfoodA=wualphaenvfood[(wualphaenvfood$cage_treatment == "Apple"),]
wualphaenvfoodB=wualphaenvfood[(wualphaenvfood$cage_treatment == "Bloom"),]
wualphaenvfoodAT=wualphaenvfood[(wualphaenvfood$cage_treatment == "Apple_AT"),]
wualphaenvfoodLB=wualphaenvfood[(wualphaenvfood$cage_treatment == "Apple_LB"),]

uualphaenvfood=uualphameta[(uualphameta$sample_type_detail == "cage food") | (uualphameta$sample_type_detail == "cage env"),]

uualphaenvfood2=uualphaenvfood[(uualphaenvfood$time_point == "T2"),]
uualphaenvfood4=uualphaenvfood[(uualphaenvfood$time_point == "T4"),]

uualphaenvfoodAB=uualphaenvfood[(uualphaenvfood$cage_treatment == "Apple") | (uualphaenvfood$cage_treatment == "Bloom"),]
uualphaenvfoodATLB=uualphaenvfood[(uualphaenvfood$cage_treatment == "Apple_AT") | (uualphaenvfood$cage_treatment == "Apple_LB"),]
uualphaenvfoodA=uualphaenvfood[(uualphaenvfood$cage_treatment == "Apple"),]
uualphaenvfoodB=uualphaenvfood[(uualphaenvfood$cage_treatment == "Bloom"),]
uualphaenvfoodAT=uualphaenvfood[(uualphaenvfood$cage_treatment == "Apple_AT"),]
uualphaenvfoodLB=uualphaenvfood[(uualphaenvfood$cage_treatment == "Apple_LB"),]

## wild flies + cage food + cage env
uualphawildcageenv=uualphameta[(uualphameta$sample_type_detail == "flies") | (uualphameta$sample_type_detail == "cage food") | (uualphameta$sample_type_detail == "cage env"),]

uualphawildcageenvAB=uualphawildcageenv[(uualphawildcageenv$cage_treatment == "Apple") | (uualphawildcageenv$cage_treatment == "Bloom"),]
uualphawildcageenvATLB=uualphawildcageenv[(uualphawildcageenv$cage_treatment == "Apple_AT") | (uualphawildcageenv$cage_treatment == "Apple_LB"),]
uualphawildcageenvA=uualphawildcageenv[(uualphawildcageenv$cage_treatment == "Apple"),]
uualphawildcageenvB=uualphawildcageenv[(uualphawildcageenv$cage_treatment == "Bloom"),]
uualphawildcageenvAT=uualphawildcageenv[(uualphawildcageenv$cage_treatment == "Apple_AT"),]
uualphawildcageenvLB=uualphawildcageenv[(uualphawildcageenv$cage_treatment == "Apple_LB"),]

uualphawildcageenv2=uualphawildcageenv[(uualphawildcageenv$time_point == "T2"),]
uualphawildcageenv4=uualphawildcageenv[(uualphawildcageenv$time_point == "T4"),]

wualphawildcageenv=wualphameta[(wualphameta$sample_type_detail == "flies") | (wualphameta$sample_type_detail == "cage food") | (wualphameta$sample_type_detail == "cage env"),]

wualphawildcageenv2=wualphawildcageenv[(wualphawildcageenv$time_point == "T2"),]
wualphawildcageenv4=wualphawildcageenv[(wualphawildcageenv$time_point == "T4"),]

wualphawildcageenvAB=wualphawildcageenv[(wualphawildcageenv$cage_treatment == "Apple") | (wualphawildcageenv$cage_treatment == "Bloom"),]
wualphawildcageenvATLB=wualphawildcageenv[(wualphawildcageenv$cage_treatment == "Apple_AT") | (wualphawildcageenv$cage_treatment == "Apple_LB"),]
wualphawildcageenvA=wualphawildcageenv[(wualphawildcageenv$cage_treatment == "Apple"),]
wualphawildcageenvB=wualphawildcageenv[(wualphawildcageenv$cage_treatment == "Bloom"),]
wualphawildcageenvAT=wualphawildcageenv[(wualphawildcageenv$cage_treatment == "Apple_AT"),]
wualphawildcageenvLB=wualphawildcageenv[(wualphawildcageenv$cage_treatment == "Apple_LB"),]

## wild flies + f1s + f2s 
uualphawildF1sF2s=uualphameta[(uualphameta$sample_type_detail == "flies") | (uualphameta$sample_type_detail == "dry stored f1s") | (uualphameta$sample_type_detail == "dry stored f2s"),]

uualphawildF1sF2sAB=uualphawildF1sF2s[(uualphawildF1sF2s$cage_treatment == "Apple") | (uualphawildF1sF2s$cage_treatment == "Bloom"),]
uualphawildF1sF2sATLB=uualphawildF1sF2s[(uualphawildF1sF2s$cage_treatment == "Apple_AT") | (uualphawildF1sF2s$cage_treatment == "Apple_LB"),]
uualphawildF1sF2sA=uualphawildF1sF2s[(uualphawildF1sF2s$cage_treatment == "Apple"),]
uualphawildF1sF2sB=uualphawildF1sF2s[(uualphawildF1sF2s$cage_treatment == "Bloom"),]
uualphawildF1sF2sAT=uualphawildF1sF2s[(uualphawildF1sF2s$cage_treatment == "Apple_AT"),]
uualphawildF1sF2sLB=uualphawildF1sF2s[(uualphawildF1sF2s$cage_treatment == "Apple_LB"),]

uualphawildF1sF2s2=uualphawildF1sF2s[(uualphawildF1sF2s$time_point == "T2"),]
uualphawildF1sF2s4=uualphawildF1sF2s[(uualphawildF1sF2s$time_point == "T4"),]

wualphawildF1sF2s=wualphameta[(wualphameta$sample_type_detail == "flies") | (wualphameta$sample_type_detail == "dry stored f1s") | (wualphameta$sample_type_detail == "dry stored f2s"),]

wualphawildF1sF2s2=wualphawildF1sF2s[(wualphawildF1sF2s$time_point == "T2"),]
wualphawildF1sF2s4=wualphawildF1sF2s[(wualphawildF1sF2s$time_point == "T4"),]

wualphawildF1sF2sAB=wualphawildF1sF2s[(wualphawildF1sF2s$cage_treatment == "Apple") | (wualphawildF1sF2s$cage_treatment == "Bloom"),]
wualphawildF1sF2sATLB=wualphawildF1sF2s[(wualphawildF1sF2s$cage_treatment == "Apple_AT") | (wualphawildF1sF2s$cage_treatment == "Apple_LB"),]
wualphawildF1sF2sA=wualphawildF1sF2s[(wualphawildF1sF2s$cage_treatment == "Apple"),]
wualphawildF1sF2sB=wualphawildF1sF2s[(wualphawildF1sF2s$cage_treatment == "Bloom"),]
wualphawildF1sF2sAT=wualphawildF1sF2s[(wualphawildF1sF2s$cage_treatment == "Apple_AT"),]
wualphawildF1sF2sLB=wualphawildF1sF2s[(wualphawildF1sF2s$cage_treatment == "Apple_LB"),]

##cage food + wild flies + f1s
uualphawildcagef1s=uualphameta[(uualphameta$sample_type_detail == "flies") | (uualphameta$sample_type_detail == "cage food") | (uualphameta$sample_type_detail == "dry stored f1s"),]

uualphawildcagef1sAB=uualphawildcagef1s[(uualphawildcagef1s$scage_treatment == "Apple") | (uualphawildcagef1s$cage_treatment == "Bloom"),]
uualphawildcagef1sATLB=uualphawildcagef1s[(uualphawildcagef1s$cage_treatment == "Apple_AT") | (uualphawildcagef1s$cage_treatment == "Apple_LB"),]
uualphawildcagef1sA=uualphawildcagef1s[(uualphawildcagef1s$cage_treatment == "Apple"),]
uualphawildcagef1sB=uualphawildcagef1s[(uualphawildcagef1s$cage_treatment == "Bloom"),]
uualphawildcagef1sAT=uualphawildcagef1s[(uualphawildcagef1s$cage_treatment == "Apple_AT"),]
uualphawildcagef1sLB=uualphawildcagef1s[(uualphawildcagef1s$cage_treatment == "Apple_LB"),]

uualphawildcagef1s2=uualphawildcagef1s[(uualphawildcagef1s$time_point == "T2"),]
uualphawildcagef1s4=uualphawildcagef1s[(uualphawildcagef1s$time_point == "T4"),]

wualphawildcagef1s=wualphameta[(wualphameta$sample_type_detail == "flies") | (wualphameta$sample_type_detail == "cage food") | (wualphameta$sample_type_detail == "dry stored f1s"),]

wualphawildcagef1sAB=wualphawildcagef1s[(wualphawildcagef1s$cage_treatment == "Apple") | (wualphawildcagef1s$cage_treatment == "Bloom"),]
wualphawildcagef1sATLB=wualphawildcagef1s[(wualphawildcagef1s$cage_treatment == "Apple_AT") | (wualphawildcagef1s$cage_treatment == "Apple_LB"),]
wualphawildcagef1sA=wualphawildcagef1s[(wualphawildcagef1s$cage_treatment == "Apple"),]
wualphawildcagef1sB=wualphawildcagef1s[(wualphawildcagef1s$cage_treatment == "Bloom"),]
wualphawildcagef1sAT=wualphawildcagef1s[(wualphawildcagef1s$cage_treatment == "Apple_AT"),]
wualphawildcagef1sLB=wualphawildcagef1s[(wualphawildcagef1s$cage_treatment == "Apple_LB"),]

wualphawildcagef1s2=wualphawildcagef1s[(wualphawildcagef1s$time_point == "T2"),]
wualphawildcagef1s4=wualphawildcagef1s[(wualphawildcagef1s$time_point == "T4"),]


###VIZ PCAS

### with more than one cage treatment 
wualphawildF1sAB%>%
  filter (!cage_treatment == 'Bloom') %>%
  ggplot(aes(x=PC1, y=PC3, color= cage_treatment, shape= sample_type_detail)) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(19,1,17),name="Sample Type") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Cage Treatment")
#### THIS PCA MA ####
uualphaenvfood%>%
  filter (!cage_treatment == 'Bloom') %>% 
  filter (!cage_treatment == 'base Bloom') %>% 
  filter (!sample_type_detail == 'cage env') %>% 
  filter (!time_point == 'T5') %>%   
  ggplot(aes(x=PC1, y=PC2, color= cage_treatment, shape=interaction(time_point, sample_type_detail))) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(19,8),name="Sample Type", labels = c("T2","T4")) + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Apple_LB", "Apple_AT"),
                      labels = c("Control Populations","LB+ Populations", "AT+ Populations"),
                      values = c("#848FA2","#058ED9","#CC2D35"))+
  theme_bw()  

#### THIS PCA AB ####
View(uualphaF1s)
a=wualphawild%>%
  filter (!cage_treatment == 'Apple_AT') %>% 
  filter (!cage_treatment == 'Apple_LB') %>% 
  filter (!cage_treatment == 'base Bloom') %>% 
  filter (!sample_type_detail == 'cage env') %>% 
  filter (!time_point == 'T5') %>% 
  filter (!time_point == 'T0') %>% 
  ggplot(aes(x=PC1, y=PC2, color= cage_treatment, shape=interaction(time_point, sample_type_detail))) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + 
  stat_ellipse(type='t',size =1)+#alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  #facet_wrap(. ~sample_type_detail)+
  scale_shape_manual(values=c(19,8),name="Sample Type", labels = c("T2","T4")) + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple Populations","Bloomington Populations"),
                      values = c("#CC2D35","black"))+
  theme_bw()  

a

#View(uualphawildcageLB)

#### THIS PCA MA

uualphaenvfood%>%
  filter (!cage_treatment == 'Apple_LB') %>% 
  filter (!cage_treatment == 'Apple_AT') %>% 
  #filter (!cage_treatment == 'base Bloom') %>% 
  filter (!sample_type_detail == 'cage env') %>% 
  filter (!time_point == 'T5') %>%   
  ggplot(aes(x=PC1, y=PC2, color= cage_treatment, shape=interaction(time_point, sample_type_detail))) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(19,8),name="Sample Type", labels = c("T2","T4")) + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Apple","Bloom"),
                      values = c("red","black"))+
  theme_bw()  

View(uualphawildcageLB)

wualphawildF1sAb%>%
  ggplot(aes(x=PC1, y=PC3, color= time_point, shape= sample_type_detail)) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(7,17),name="Sample Type") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Time Point")


###only one cage treatment
View (wualphaF1sB)

wualphawild%>%
  ggplot(aes(x=PC1, y=PC2, color= time_point, shape=interaction( sample_type_detail))) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_discrete(name="Sample Type") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Time Point")





####Individual PCAs (less data points)###################################################











####heatmaps#############################1426######################
library(tidyverse)
library(qiime2R)
library(ggpubr)

#metadata<-read_q2metadata("orchard2016smetadata.csv")
SVs<-read_qza("table-no-wolbachia-exact.qza")$data
SVs<-apply(SVs, 2, function(x) x/sum(x)*100) #convert to percent
taxonomy<-read_qza("correcttaxonomy.qza")$data %>% parse_taxonomy()
View(taxonomy)
##long to run 
taxasumgenus<-summarize_taxa(SVs, taxonomy)$Genus
taxasumfamily<-summarize_taxa(SVs, taxonomy)$Family
taxasumorder<-summarize_taxa(SVs, taxonomy)$Order
taxasumclass<-summarize_taxa(SVs, taxonomy)$Class
taxasumking<-summarize_taxa(SVs, taxonomy)$Kingdom
View(alphadata)
#subsetting alphadata
View(alphadata)
alphadata$treatment_timepoint=paste(alphadata$cage_treatment, alphadata$time_point)
alphadatanc = subset(alphadata, alphadata$control == "no" , na.rm=TRUE)
founder=subset(alphadatanc,alphadatanc$time_point== "T0", na.rm=TRUE)
wild=subset(alphadata, alphadata$sample_type_detail == "flies", na.rm=TRUE)
wild2=subset(wild, wild$time_point== "T2", na.rm=TRUE)
wild4=subset(wild, wild$time_point== "T4", na.rm=TRUE)
food=subset(alphadata, alphadata$sample_type_detail == "cage food")
env=subset(alphadata, alphadata$sample_type_detail == "cage env")
f1s=subset(alphadata, alphadata$sample_type_detail == "dry stored f1s")
f2s=subset(alphadata, alphadata$sample_type_detail == "dry stored f2s")
nobloom=subset(alphadata, alphadata$cage_treatment%in% c("Apple", "Apple_AT", "Apple_LB"))
View(f2s24)
alphadata24=subset(nobloom, nobloom$time_point %in% c("T2", "T4"))
alphadata2and4=subset(alphadata, alphadata$time_point %in% c("T2", "T4"))
wild24=subset(alphadata24, alphadata24$sample_type_detail == "flies")
food24=subset(alphadata24, alphadata24$sample_type_detail == "cage food")
env24=subset(alphadata24, alphadata24$sample_type_detail == "cage env")
f1s24=subset(alphadata24, alphadata24$sample_type_detail == "dry stored f1s")
f2s24=subset(alphadata2and4, alphadata2and4$sample_type_detail == "dry stored f2s")
apple24=subset(alphadata24, alphadata24$cage_treatment == "Apple")
applewild24=subset(apple24, apple24$sample_type_detail %in% c("flies"))
applefood24=subset(apple24, apple24$sample_type_detail %in% c("cage food"))
View(apple24)

View(founder)
View(taxasumgenus)
####Wild figures####

taxasumgenus %>%
  as.data.frame(t(taxasumgenus)) %>%
  filter[founder$SampleID, ] %>%
  as.data.frame(t(taxasumgenus))
  
##Genus##
##all founder
t.taxasumgenus=as.data.frame(t(taxasumgenus))
View(t.taxasumgenus)
t.foundertaxasumgenus=t.taxasumgenus[founder$SampleID, ]
View(t.foundertaxasumgenus)
foundertaxasumgenus=as.data.frame(t(t.foundertaxasumgenus))

founderall=taxa_heatmap(foundertaxasumgenus, founder, "study_group", normalize = TRUE, ntoplot = 10)
founderall

##all wild
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.wildtaxasumgenus=t.taxasumgenus[wild24$SampleID, ]
wildtaxasumgenus=as.data.frame(t(t.wildtaxasumgenus))

##allf1s
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.F1staxasumgenus=t.taxasumgenus[f1s24$SampleID, ]
F1staxasumgenus=as.data.frame(t(t.F1staxasumgenus))

##allf2s
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.F2staxasumgenus=t.taxasumgenus[f2s24$SampleID, ]
F2staxasumgenus=as.data.frame(t(t.F2staxasumgenus))

genuswild24all=taxa_heatmap(wildtaxasumgenus, wild24, "treatment_timepoint", normalize = TRUE, ntoplot = 15)
genuswild24all

genusF124all=taxa_heatmap(F1staxasumgenus, f1s24, "treatment_timepoint", normalize = TRUE, ntoplot = 15)
genusF124all

genusF224all=taxa_heatmap(F2staxasumgenus, f2s24, "treatment_timepoint", normalize = TRUE, ntoplot = 15)
genusF224all

#applewild


t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.AWtaxasumgenus=t.taxasumgenus[applewild24$SampleID, ]
AWtaxasumgenus=as.data.frame(t(t.AWtaxasumgenus))

View(applewild24)
genusAPPLEwild=taxa_heatmap(AWtaxasumgenus, applewild24, "time_point", normalize = TRUE, ntoplot = 15)
genusAPPLEwild


#applefood
t.taxasumgenus=as.data.frame(t(taxasumgenus))
t.AFtaxasumgenus=t.taxasumgenus[applefood24$SampleID, ]
AFtaxasumgenus=as.data.frame(t(t.AFtaxasumgenus))

View(applefood24)
genusAPPLEfood=taxa_heatmap(AFtaxasumgenus, applefood24, "time_point", normalize = TRUE, ntoplot = 15)
genusAPPLEfood

library(cowplot)
plot_grid(genusAPPLEwild,genusAPPLEfood, labels = c('Wild Flies', 'Cage Food'), ncol = 1)

##Family##

##all wild
t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.wildtaxasumfamily=t.taxasumfamily[wild24$SampleID, ]
wildtaxasumfamily=as.data.frame(t(t.wildtaxasumfamily))

##all F1s
t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.f1taxasumfamily=t.taxasumfamily[f1s24$SampleID, ]
f1taxasumfamily=as.data.frame(t(t.f1taxasumfamily))

##all f2s
t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.f2taxasumfamily=t.taxasumfamily[f2s24$SampleID, ]
f2taxasumfamily=as.data.frame(t(t.f2taxasumfamily))

View(f2s24)


familywild24all=taxa_heatmap(wildtaxasumfamily, wild24, "treatment_timepoint", ntoplot = 10)
familywild24all

familyf124all=taxa_heatmap(f1taxasumfamily, f1s24, "treatment_timepoint", normalize = TRUE,  ntoplot = 10)
familyf124all

familyf224all=taxa_heatmap(f2taxasumfamily, f2s24, "treatment_timepoint", normalize = TRUE,  ntoplot = 10)
familyf224all

#applewildfam
t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.AWtaxasumfamily=t.taxasumfamily[applewild24$SampleID, ]
AWtaxasumfamily=as.data.frame(t(t.AWtaxasumfamily))

View(applewild24)
familyAPPLEwild=taxa_heatmap(AWtaxasumfamily, applewild24, "time_point", normalize = TRUE, ntoplot = 15)
familyAPPLEwild


#applefoodfam
t.taxasumfamily=as.data.frame(t(taxasumfamily))
t.AFtaxasumfamily=t.taxasumfamily[applefood24$SampleID, ]
AFtaxasumfamily=as.data.frame(t(t.AFtaxasumfamily))

View(applefood24)
familyAPPLEfood=taxa_heatmap(AFtaxasumfamily, applefood24, "time_point", normalize = TRUE, ntoplot = 15)
familyAPPLEfood

library(cowplot)
plot_grid(familyAPPLEwild,familyAPPLEfood, labels = c('Wild Flies', 'Cage Food'), ncol = 1)

##Order##

##all wild
t.taxasumorder=as.data.frame(t(taxasumorder))
t.wildtaxasumorder=t.taxasumorder[wild24$SampleID, ]
wildtaxasumorder=as.data.frame(t(t.wildtaxasumorder))

##all F1s
t.taxasumorder=as.data.frame(t(taxasumorder))
t.f1taxasumorder=t.taxasumorder[f1s24$SampleID, ]
f1taxasumorder=as.data.frame(t(t.f1taxasumorder))


orderwild24all=taxa_heatmap(wildtaxasumorder, wild24, "treatment_timepoint", ntoplot = 10)
orderwild24all

orderf124all=taxa_heatmap(f1taxasumorder, f1s24, "treatment_timepoint", normalize = TRUE,  ntoplot = 10)
orderf124all

##class##

##all wild
t.taxasumclass=as.data.frame(t(taxasumclass))
t.wildtaxasumclass=t.taxasumclass[wild24$SampleID, ]
wildtaxasumclass=as.data.frame(t(t.wildtaxasumclass))

##all F1s
t.taxasumclass=as.data.frame(t(taxasumclass))
t.f1taxasumclass=t.taxasumclass[f1s24$SampleID, ]
f1taxasumclass=as.data.frame(t(t.f1taxasumclass))


classwild24all=taxa_heatmap(wildtaxasumclass, wild24, "treatment_timepoint", ntoplot = 10)
classwild24all

classf124all=taxa_heatmap(f1taxasumclass, f1s24, "treatment_timepoint", normalize = TRUE,  ntoplot = 10)
classf124all

##kingdom##

##all wild
t.taxasumking=as.data.frame(t(taxasumking))
t.wildtaxasumking=t.taxasumking[wild24$SampleID, ]
wildtaxasumking=as.data.frame(t(t.wildtaxasumking))
View(wildtaxasumking)
##all F1s
t.taxasumking=as.data.frame(t(taxasumking))
t.f1taxasumking=t.taxasumking[f1s24$SampleID, ]
f1taxasumking=as.data.frame(t(t.f1taxasumking))


kingwild24all=taxa_heatmap(wildtaxasumking, wild24, "treatment_timepoint", ntoplot = 10)
kingwild24all

classf124all=taxa_heatmap(f1taxasumclass, f1s24, "treatment_timepoint", normalize = TRUE,  ntoplot = 10)
classf124all


##all F2s
t.taxasumking=as.data.frame(t(taxasumking))
t.f2taxasumking=t.taxasumking[f2s24$SampleID, ]
f2taxasumking=as.data.frame(t(t.f1taxasumking))


kingwild24all=taxa_heatmap(wildtaxasumking, wild24, "treatment_timepoint", ntoplot = 10)
kingwild24all

classf124all=taxa_heatmap(f1taxasumclass, f1s24, "treatment_timepoint", normalize = TRUE,  ntoplot = 10)
classf124all

######ALTERNATIVE APPROACH top 15 read abundance######
library(tidyverse)
library(dplyr)
library(qiime2R)

metadata<-alphadata
alphadataSVs<-read_qza("table-no-wolbachia-exact.qza")$data
taxonomy<-read_qza("correcttaxonomy.qza")$data
library(tibble)
taxonomy <- tibble::rownames_to_column(taxonomy, "Feature.ID")

SVs<-apply(SVs, 2, function(x) x/sum(x)*100) #convert to percent

SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(15, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table
View(SVsToPlot)
View(taxonomy)
View(SVs)
View(tim)
View(metadata)
BIGSVs=SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  #summarize(Abundance=sum(Abundance)) %>%
  left_join(alphadata, by= "SampleID") %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy, by= "Feature.ID") %>%
  mutate(Feature=paste(Feature.ID, Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature))  # trim out leading text from taxonomy string
View(BIGSVs)


SmallSVs=BIGSVs[!(BIGSVs$Feature.ID == "Remainder"),]
test=SmallSVs[sample(nrow(SmallSVs), 10), ]
test$Feature=gsub(".*Bacteria;","",as.character(test$Feature)) ####fixes trimming issue of feature variable

SmallSVs=BIGSVs[!(BIGSVs$Feature.ID == "Remainder"),]##remove remainder
SmallSVs$Feature=gsub(".*Bacteria;","",as.character(SmallSVs$Feature))

BIGSVs$Feature=gsub(".*Bacteria;","",as.character(BIGSVs$Feature))


####Subsetting new BIGSVs####
wildSVs=subset(BIGSVs, BIGSVs$sample_type_detail == "flies", na.rm=TRUE)
wild2SVs=subset(wildSVs, wildSVs$time_point== "T2", na.rm=TRUE)
wild4SVs=subset(wildSVs, wildSVs$time_point== "T4", na.rm=TRUE)
all24SVs=subset(BIGSVs, BIGSVs$time_point %in% c("T2", "T4"), na.rm=TRUE)
apple24SVs=subset(all24SVs,all24SVs$cage_treatment == "Apple", na.rm=TRUE)
applewildfood24SVs=subset(apple24SVs, apple24SVs$sample_type_detail %in% c("flies", "cage food"))
View(applewildfood24SVs)

foodSVs=subset(BIGSVs, BIGSVs$sample_type_detail == "cage food")
food2SVs=subset(foodSVs, foodSVs$time_point== "T2", na.rm=TRUE)
food4SVs=subset(foodSVs, foodSVs$time_point== "T4", na.rm=TRUE)

envSVs=subset(BIGSVs, BIGSVs$sample_type_detail == "cage env")
founderSVs=subset(BIGSVs, BIGSVs$time_point == "T0")

f1sSVs=subset(BIGSVs, BIGSVs$sample_type_detail == "dry stored f1s")
f1sSVsAPATLB=subset(f1sSVs, f1sSVs$cage_treatment %in% c("Apple", "Apple_AT", "Apple_LB"))
f1sSVsAPATLB24=subset(f1sSVsAPATLB, f1sSVsAPATLB$time_point %in% c("T2", "T4"))


f2sSVst=subset(BIGSVs, BIGSVs$sample_type_detail == "dry stored f2s")
f2sSVs=subset(f2sSVst, f2sSVst$cage_treatment %in% c("Apple_AX", "Axenic_AT", "Axenic_LB"))
alphawildF1s <-alphadata[(alphadata$sample_type_detail == "flies") | (alphadata$sample_type_detail == "dry stored f1s"),]

##Apple T2 /T4 alone##
AppleT2SVs<-SmallSVs[(SmallSVs$time_point == "T2") & (SmallSVs$cage_treatment == "Apple"),]
AppleT2SVs=AppleT2SVs[complete.cases(AppleT2SVs), ]
AppleT4SVs<-SmallSVs[(SmallSVs$time_point == "T4") & (SmallSVs$cage_treatment == "Apple"),]
AppleT4SVs=AppleT4SVs[complete.cases(AppleT4SVs), ]

##AT T2 /T4 alone##
AppleATT2SVs<-SmallSVs[(SmallSVs$time_point == "T2") & (SmallSVs$cage_treatment %in% c("Apple_AT", "Axeinc_AT")),]
AppleATT2SVs=AppleATT2SVs[complete.cases(AppleATT2SVs), ]
AppleATT4SVs<-SmallSVs[(SmallSVs$time_point == "T4") & (SmallSVs$cage_treatment %in% c("Apple_AT", "Axeinc_AT")),]
AppleATT4SVs=AppleATT4SVs[complete.cases(AppleATT4SVs), ]

##LB T2 /T4 alone##
AppleLBT2SVs<-SmallSVs[(SmallSVs$time_point == "T2") & (SmallSVs$cage_treatment %in% c("Apple_LB", "Axeinc_LB")),]
AppleLBT2SVs=AppleLBT2SVs[complete.cases(AppleLBT2SVs), ]
AppleLBT4SVs<-SmallSVs[(SmallSVs$time_point == "T4") & (SmallSVs$cage_treatment %in% c("Apple_LB", "Axeinc_LB")),]
AppleLBT4SVs=AppleLBT4SVs[complete.cases(AppleLBT4SVs), ]

##Bloom T2 /T4 alone##
BloomT2SVs<-SmallSVs[(SmallSVs$time_point == "T2") & (SmallSVs$cage_treatment == "Bloom"),]
BloomT2SVs=BloomT2SVs[complete.cases(BloomT2SVs), ]
BloomT4SVs<-SmallSVs[(SmallSVs$time_point == "T4") & (SmallSVs$cage_treatment == "Bloom"),]
BloomT4SVs=BloomT4SVs[complete.cases(BloomT4SVs), ]

##wild flies treatmentxtimepoint## 

a=ggplot(wildSVs, aes(x=study_group, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
  theme(legend.position = "bottom") +
  ggtitle("Wild") +
  xlab("Cage Treatment x Season")+ 
  ylab("15 Most Abundant Bacterial Reads")+
  scale_fill_viridis_c(option = "plasma" , name="Abundance",  )

a 
 #ggsave("wildflyheatmap.pdf  ", height=4, width=10, device="pdf") # save a PDF ? inches by ?  inches



##cage food treatmentxtimepoint## 

  ggplot(foodSVs, aes(x=study_group, y=Feature, fill=NormAbundance)) +
    geom_tile() +
    facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
    theme_q2r() +
    theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
    theme(legend.position = "bottom") +
    ggtitle("Cage Food Samples") +
    xlab("Cage Treatment x Season")+ 
    ylab("15 Most Abundant Bacterial Reads")+
    scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
 
   ggsave("cagefoodheatmap.pdf  ", height=4, width=10, device="pdf") # save a PDF ? inches by ?  inches


##cage env treatmentxtimepoint## 
   
   
   ggplot(envSVs, aes(x=study_group, y=Feature, fill=NormAbundance)) +
     geom_tile() +
     facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("Cage ENV Samples") +
     xlab("Cage Treatment x Season")+ 
     ylab("15 Most Abundant Bacterial Reads")+
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("cageenvheatmap.pdf  ", height=4, width=7, device="pdf") # save a PDF ? inches by ?  inches
   
   
##f1s treatment x timepoint###
f1sSVs=subset(f1sSVs, f1sSVs$time_point %in% c("T2","T4"),)
   ggplot(f1sSVs, aes(x=study_group, y=Feature, fill=NormAbundance)) +
     geom_tile() +
     facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("F1 Samples") +
     xlab("Cage Treatment x Season") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("F1heatmap24.pdf  ", height=4, width=10, device="pdf") # save a PDF ? inches by ?  inches
   
   
##f1s treatment x timepoint###
   
   ggplot(f2sSVs, aes(x=study_group, y=Feature, fill=NormAbundance)) +
     geom_tile() +
     facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("F2 Samples") +
     xlab("Cage Treatment x Season") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("F2heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   ##all summer apple###

   ggplot(AppleT2SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("Apple T2 Samples") +
     xlab("Cage Treatment x Season") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("appleT2heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   
   ##all fall apple###
   
   ggplot(AppleT4SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("Apple T4 Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("appleT4heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   
   ##all summer apple AT ###
   
   ggplot(AppleATT2SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("AT T2 Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("appleATT2heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   
   ##all fall AT ###
   
   ggplot(AppleATT4SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("AT T4 Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("appleATT4heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   ##all summer apple LB###
   
   ggplot(AppleLBT2SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("LB T2 Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("appleLBT2heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   
   ##all fall LB ###
   
   ggplot(AppleLBT4SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("LB T4 Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("appleLBT4heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   ##all summer bloom###
   
   ggplot(BloomT2SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("Bloom T2 Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("BloomT2heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
   
   ##all fall bloom ###
   
   ggplot(BloomT4SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`sample_type_detail`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("BloomT4 Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("BLoomT4heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
  
   ggplot(founderSVs, aes(x=sample_type_detail, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`study_group`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("Founder Samples") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggsave("BLoomT4heatmap.pdf  ", height=4, width=8, device="pdf") # save a PDF ? inches by ?  inches
   
####apple wild food####
   ggplot(applewildfood24SVs, aes(x=study_group, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~interaction(time_point,sample_type_detail), scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("Apple wild food ") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
   ggplot(applewildfood24SVs, aes(x=study, y=Feature, fill=NormAbundance), na.rm=TRUE) +
     geom_tile() +
     facet_grid(~`study_group`, scales="free_x") +
     theme_q2r() +
     theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
     theme(legend.position = "bottom") +
     ggtitle("Apple wild food ") +
     xlab("Sample Type") + 
     ylab("15 Most Abundant Bacterial Reads") +
     scale_fill_viridis_c(option = "plasma" , name="Abundance",  )
   
####volcano#### 
   ##need to run q2-aldex diff abundance analysis in qiime
   library(tidyverse)
   library(qiime2R)
   library(ggrepel) # for offset labels
   devtools::install_github("GuangchuangYu/ggtree")
   library(ggtree) # for visualizing phylogenetic trees
   library(ape) # for manipulating phylogenetic trees
   
   metadata<-alphadata
   SVs<-read_qza("table-no-wolbachia-exact.qza")$data
   results<-read_qza("differentials.qza")$data #### dont have this
   taxonomy<-read_qza("correcttaxonomy.qza")$data
   tree<-read_qza("rooted-tree.qza")$data
   ##Volcano Plot
   results %>%
     left_join(taxonomy) %>%
     mutate(Significant=if_else(we.eBH<0.1,TRUE, FALSE)) %>%
     mutate(Taxon=as.character(Taxon)) %>%
     mutate(TaxonToPrint=if_else(we.eBH<0.1, Taxon, "")) %>% #only provide a label to signifcant results
     ggplot(aes(x=diff.btw, y=-log10(we.ep), color=Significant, label=TaxonToPrint)) +
     geom_text_repel(size=1, nudge_y=0.05) +
     geom_point(alpha=0.6, shape=16) +
     theme_q2r() +
     xlab("log2(fold change)") +
     ylab("-log10(P-value)") +
     theme(legend.position="none") +
     scale_color_manual(values=c("black","red"))
   ggsave("volcano.pdf", height=3, width=3, device="pdf")
   
####per-feature abundance####
   View(SVs)
   clr<-apply(log2(SVs+0.5), 2, function(x) x-mean(x))
   
   View(alphadata24)
   View(taxonomy)
   View(clr)
   
   LacBrevFID= c("6dc086bb99734438fbeb2324e96c9f51",
                 "0ae6079001a4db920e51b95470b38f4b",
                 "357b9856061a7ca06c2c4dcc212f5cc1",
                 "4f7170de2752a62232a1923f04c046a4",
                 "c2a2284f9cce143499ab156ebb28e3ce",
                 "ceb32fddb3b21ebd5d3ce86297bd951f",
                 "efcc1b159b1927cdb218803ea4a4b22c",
                 "06e098f2352f50dee72776a6371f884c",
                 "2666688937208db0d0aef1c1185125a4",
                 "6ab936d69229e7655c60144c3cb52469",
                 "a44090ab73b9ecc49f3badc35d60ccb5",
                 "70282b3ede979f80398f175e6c89869d",
                 "aacda0d1ea0473bd03152a7b1dcacab8",
                 "c0ced4099187985720aed7630b4bdd32",
                 "06794eab8b662c5add16c12f45ec5ee9",
                 "fc042e344b3bccffc927f9a2d376d82c",
                 "c3d072802c8a5bc2ef1c89c58333a7e3",
                 "7050eb9051a55720145a040c069da1d5",
                 "86bc687f9a33c723e0682f13dee4a78e")
   
   AcetoFID=c("4d0797327ef46bf483334d7bfa376255",
              "f2da1736c3021585914c1fff306f3f2e",
              "6459a3b47a2ea884885b3095088480f1",
              "c882cb61faf20ec04253b3443b6f72bd",
              "c9d2b6e9be51e5fd6a62c3885c8b5f42",
              "b3a92e475daf848eb76f35a531306a1b",
              "8172bd0a1889f5538a4f901fa2be210f",
              "0b6b0090a979decfdeb4f727d0c06252",
              "e0b6aab721feda59dea0bdb10cb56c72",
              "5b4fac426fb13a5d4f7ad8890c578fb7",
              "e685172159e7edbe8ee1d928328994a9",
              "cd792bae2aaa747a26d1cf2fb2545289",
              "099c6bf73c4083e1164ff2c7977725ca",
              "7a2853a64b6d1995f970ddb3a17670e4",
              "6d95f1b6311a87772ba6daf4ddd61ea2",
              "bb6b9883c1db30b019ea457f1a250bdf",
              "5628f17221f8817ca3b4462cc19d6070")
   
  GluconFID=c("6c0be6da4b72648fd28d31272603c9eb",
              "074dfcbdcc5db6da81d9a3da3895a85e",
              "93f7ff7537231b1ae7a80ccc1848f42f",
              "d6bb0a532019f7d722debcd2189b8c8b",
              "deb86ea427ab80f68c3dfcb354544546",
              "a44e184808af4ff62921951f74a20c77",
              "c23a349b8c516ce0eab2337f4c4d7511",
              "5c30ab3776eef304c1c978bb65486175",
              "3320125f47e97864c6b76777ffe8d38d",
              "2a3a8185ceb551ad8941b5cdd012009a",
              "44096639a24f1a315b06822a33ee4cf0",
              "d4b8211a2e41f1fd13578891002b4760",
              "ab296c42ad90ea27b863b3e374e36807",
              "8ee154513cbcd7b2a713819d64227c43")   
   
  GluconacFID=c("38e8b158e5e76ffb60a34672886b1fae",
                "8e0d90a7aafa7f322e1c370659e52974",
                "f26d0844a19b803f78c0b21133692247",
                "cb4f3b1689d284c92bebbbad0361b9c6",
                "d09c015d59d9276d22de803a5db02822",
                "0130b94cbe554dc6b6fff3b2ca8b5a0d",
                "0b86c846bda2e07d22c01fa23d6a3300")   

  EnterobacterFID=c("ac048f828e1f72fdd41184a35e0a3109",
                "9c678a953d519e95bd9b0cebe447c986",
                "f25b61f59833ff29e2b6445bd55ffea1",
                "9bc0642346e4428f14dbfeae4e973dc2",
                "57c85d581be762cf58c8277f63d860e3")   
  
  LeuconoFID=c("9c692e302675ffd1d66e0039675fd533",
               "8b73486ee4e7877c6fb8872946e5b139",
               "da61627240af7b4694b0ad2d6d0dc4fd",
               "0886e2cf63b64ba95c17fd9015d9f7df",
               "fd0e61bbf618d175590c87add296081d",
               "0d5b4ca12419dbb912984cbf9da3864b",
               "38d6cf514430b8d9378d1597cb99c756",
               "766b8f1cbba090e828511d5e4cdce649",
               "dc89c8e33e16b0b26a64fc2dbd0d9ed2",
               "16b898d818d6c2dd04d8265bee95533d",
               "2b6e5197311bf50d755abe215c4aaa89")
  
  ProvidFID=c("89b41bc334d91550bc522af0935fbc6a",
              "be42fe2fef09357e9e400a8031f8a964",
              "852a120eb92838bbaa6ac371a64210e6",
              "21662883b293cbddfddbae6f8b68b8aa",
              "af3043ee52654eb929f19748e86a0872",
              "378fe013916fe27e06f65983d2be6515",
              "df4f395c892f7d5ee22ad80dae01b3d6",
              "0f619cac368193b836fdab9c07b7bb2a",
              "0c34e7c70ec76cb89cb44950590ed357",
              "2a44525fc76e080a21ca6f038a2a9835",
              "980c48da4b500c2ab72e1494ce8c4523",
              "f7e72ee62f3519a44e15af0d81733bec",
              "422b0f73fcdc6812af03b9d455f0f7a8",
              "ca284ef7b16a39068b2ff6c4b7cac695",
              "cc57529ef36fab9e39d36edf937cef73",
              "016629c04fdecb1ee82c3ced2782ca30",
              "913402cdc80eb0671134bcce3449dacc",
              "0bd58e468e3ea0a2b45dd378cdb19a38",
              "50abdb1d4b499140ae87ad4d871e138d",
              "a231e308f14fb724bef89d7b77098a93",
              "7646c0231eb527071de1770cbb18ad9c",
              "903a4a6c43cbda834e9c99142ffd2214",
              "ec6fa95128685ebf9eaf3955f67c6540",
              "e2c1e33cc5ba11af178ff319746f2e0e",
              "e2c1e33cc5ba11af178ff319746f2e0e",
              "e807711385f705eecbed9457fa14b390",
              "7d0c67decaaab73daa61813ca979fb32",
              "cf4821ddc573137ece371c1b9131af12",
              "310768271c7fe31e9dd410b5a2610dc7",
              "d954aadb3ae04931f2b2aaf099702701",
              "2dd53243f6c7d2f067eb127995bc44b0")
  
    RhizoFID=c("5537e5f01e2eab34cb9b8dba44dcb864",
               "3761bf52bf2a7bb8c5c830ba3bdcbb45",
               "90bf8f4ab13ddf8fcecfc8d83b73adb2",
               "582ec96accfa1991fda2b8b9ca58af4e",
               "24a3f89d9ecc21ee75157b03f2305a6b",
               "4bc8e5c42d7aa8202446bde1177e0628",
               "47259051b049975dcbc96b5c68e72b9f",
               "f40217f677f097065be4a6693a3589ea",
               "43905d7337c08e96b73512b56b751cd5",
               "5f169e7ac2e4b1c0a85f20dfa43d187e",
               "25eedba7cf428b9a4019ebb629bc2517",
               "e0eccc4a63a52c29be3f1fed968eb005",
               "e33dde9ec81e3b9ed4bd449253b41ec8",
               "4fa1ce3ce05884578ce40d8e0db4ba72",
               "97ea73eb9e1cbe84b0a7188d6aa6cc8a")
  
  
  
  View(alphadata24)
  
   clr %>%
     as.data.frame() %>%
     rownames_to_column("Feature.ID") %>%
     gather(-Feature.ID, key=SampleID, value=CLR) %>%
     filter(Feature.ID %in% RhizoFID) %>%
     left_join(alphadata24) %>%
     filter(`sample_type_detail`=="cage food") %>%
     ggplot(aes(x=cage_treatment, y=CLR, fill=time_point)) +
     stat_summary(geom="bar", color="black", position = "dodge") +
     #geom_jitter(width=0.2, height=0, shape=21) +
     theme_q2r() 
   
   clr %>%
     as.data.frame() %>%
     rownames_to_column("Feature.ID") %>%
     gather(-Feature.ID, key=SampleID, value=CLR) %>%
     filter(Feature.ID %in% ProvidFID) %>%
     left_join(alphadata24) %>%
     filter(`sample_type_detail`=="dry stored f1s") %>%
     ggplot(aes(x=cage_treatment, y=CLR, fill=time_point)) +
     stat_summary(geom="bar", color="black", position = "dodge") +
     #geom_jitter(width=0.2, height=0, shape=21) +
     theme_q2r() 

   
   ggsave("aldexbar.pdf", height=2, width=1.5, device="pdf") 

   
   
   
#########################BACTERIAL ABUDNACE ON FOOD##############
 ##  EXPORT!!!!
View(otu.taxon_nw)
View(metadata)
View(allmeta)
View(alphadata)
view(alphafood24)### cage food wrong 
alphafood24=subset(metadata, metadata$SampleType== "Environmental control" )

a=as.data.frame(colSums(OTUTABLE_nw[,-1]))
colnames(a)[1] = "readcounts"
library(tibble)
readcounts <- tibble::rownames_to_column(a, "SampleID")

View(readcounts)

View(shannon_nw)
alphadata=alphadata %>% 
  merge(readcounts, all.y=TRUE) 

#alphadata = metadata
alphadata <- transform(alphadata, readcountsadjusted = alphadata$readcounts / alphadata$final_library_conc_ng_ul)
alphadata = subset(alphadata, !(study_day %in% c("None", NA)))

alphadatanoout <- subset(alphadata, alphadata$readcountsadjusted < 2000)
View(alphadatanoout)

alphadatanoout %>%
  as.data.frame() %>%
  filter(`sample_type_detail`== "cage food")  %>%
  ggplot() +
  #geom_boxplot(aes(x=time_point, y=observed_features, fill=cage_treatment)) +
  geom_point(aes(x=time_point , y=observed_features, group=study_group, fill=cage_treatment))

### THIS PLOT 
View(alphadatanoout)
alphadatanoout %>%
  as.data.frame() %>%
  filter(`sample_type_detail`== "flies")  %>%
  filter (!cage_treatment == 'Apple_LB') %>% 
  filter (!cage_treatment == 'Apple_AT') %>% 
  filter (!cage_treatment == 'base Bloom') %>% 
  filter (!sample_type_detail == 'cage env') %>% 
  filter (!time_point == 'T5') %>%   
  filter (!time_point == 'T0') %>%   
  ggplot(aes(x=interaction(time_point, cage_treatment), y=readcountsadjusted, fill=cage_treatment)) +
  geom_boxplot()  + 
  ylim(0,1000) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("Low Nutrition","High Nutrition"),
                      values = c("#CC2D35", "black"))+
  theme_bw()  



