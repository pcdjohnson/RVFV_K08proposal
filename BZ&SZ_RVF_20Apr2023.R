
####
# Compile of all RVF serology data for BacZoo 
# All analyses and plots for CS overview paper - human and livestock data
# 2024 updates - based on Coxiella script
# includes env covariates
###

library(lubridate)
library(reshape)
library(binom)
library(ggplot2)
library(ggmap)
library(ggrepel)
library(lme4)
library(DHARMa)
library(jtools)
library(sjPlot)
library(ROCit)
library(gridExtra)
library(gdata)

###
# Files to import data

# setwd() - location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###
# BacZoo household data (Source file: BacZoo_HHQCleaning_28Apr2017_JH.R)
H<-read.csv("data/BacZoo Households Cleaned 2018-07-18_NoCoords.csv")

# BacZoo env variable data (Source file: BacZoo_HUS_EnvironExtrat_v1_9jul21_AJM.R)
HE<-read.csv("BacZoo_HUS_EnvironExtrat_v1_2021-07-09_NoCoords.csv")

# BacZoo household coord data cleaned... from diff location in restricted file location - includes setting classification
HC<-read.csv("/Users/mojohalliday/ownCloud/Shared/Zoonoses Tanzania Restricted/BacZoo/BacZoo_HH_Coords_08July2020.csv")

# BacZoo animal data (Source file: BacZoo_Animals_Cleaning_19Dec2016_JH.R)
A<-read.csv("BacZoo Animals Cleaned 2016-12-19.csv")

# BacZoo human data (Source file: BacZoo_Humans_Cleaning_Jo_12Jul2021.R)
P<-read.csv("data/BacZoo_IndividualPeople_Cleaned_2021-07-12.csv")

# Sample data - if needed (Source file: BacZoo_Samples_Cleaning_03Jun2020.R)
ST<-read.csv("BacZoo_SampleTracking_Clean_03Jun2020.csv")
Sam<-read.csv("BacZoo_SampleList_Clean_03Jun2020.csv")

# Animal  RVF serology 
# import from WdG compile
AS<-read.csv("data/rvf_livestock_data_2719.csv")

# Human RVF serology
HS<-read.csv("data/rvf_human_data_2719.csv")

##############################################
# Merge all.. animal data
# max dataset = all animals with one or more RVF diagnostic datapoint and complete animal and household data

# H, HE and A
# trim to key fields first
H<-H[,c("HHID","IntDate",
        "Region", "District", "WardID","setting","Village",
        "keep_cattle", "keep_sheep", "keep_goats",
        "cattleTotal_n_kept","sheepTotal_n_kept", "goatsTotal_n_kept",
        "graze","contact_otherherds")]

# keep all HE
# merge all of these into H using HHID

H2<-merge(H,HE,by="HHID",all.x=TRUE)
setdiff(H$HHID,HE$HHID) # 3 HH with missing coord and therefore env data

# also merge in the coords
H3<-merge(H2,HC,by="HHID", all.x = TRUE)

#
A<-A[,c("HHID","IndID","species", "sex", "age_dentition","bcs",
        "serviceownherd", "serviceotherherd","femrepcon", 
        "deadoffspring", "deadoffspring_n",  "datedeadoffspring", 
        "retainedplac", "retained_n", "dateretainedplac",      
        "currentmilk",
        "collectmilk", "calfsuckling", "calfsucklingage_unit", 
        "calfsucklingage_number",  
        "soldmilk_pm", "soldmilk_py", "milkfreq_unit", "milkvol_unit",            
        "milkfreq_number", "milkvol_number", 
        "ecto_fleas", "ecto_lice", "ecto_ticks", "tickscore", 
        "vaccines", "vaccines_any","vaccines_which",
        "Anthrax", "PPR", "CBPP", "FMD",                      
        "CCPP", "ECF", "LSD", "RIFTVALLEYFEVER", "vaccine_other_spec")]

# merge animals into this combined HH dataset

# 19 Jun 2022 - check of abortion outcome data in this max dataset (before merges)
summary(as.factor(A$deadoffspring))
# set DK to NA

A$deadoffspring[A$deadoffspring=="DK"]<-NA
summary(as.factor(A$deadoffspring))
# summarise by species and setting - including all animals


HA<-merge(H3,A,by="HHID")
# only 1719 of 1766 merge - 
setdiff(A$IndID,HA$IndID) # pilot animals - fine to drop.... no Household data

# RVF animal serology - cELISA
head(AS)
# rename result as rvfv
AS<-plyr::rename(AS,c(result = "rvfv_result"))

HAS<-merge(HA,AS,by.x="IndID",by.y="SampleID")
# 1689 of 1719 merge
setdiff(HA$IndID,AS$SampleID) # looks like patchy missingness... runs of sample IDs..
# TODO - CHECK THESE?

###################################################################################
###################################################################################

# rename compiled animal dataset
ARVF<-HAS

##############################################
# Merge all.. human data
# max dataset = all people with Coxiella serology and complete individual and household data

# HH data as above for animals - use H3 

P<-P[,c("HHID","IntDate","IndividualID",
        "Sex", "Age_class", "ApproxAgeDOB", "AgeClassCombined",
        "Tribe","LiveVill_Unit","LiveVill_Number","primary_occupation",
        "anyoccupation_01_rancher", "anyoccupation_02_farmer", 
        "anyoccupation_03_gardener", "anyoccupation_04_mining",
        "anyoccupation_05_military", "anyoccupation_06_military",
        "anyoccupation_06_butcher", "anyoccupation_07_wildlifeworker",
        "anyoccupation_08_milksupplier", "anyoccupation_09_outdoormanual",
        "anyoccupation_10_artisan", "anyoccupation_11_veterinarian",
        "anyoccupation_12_housewife", "anyoccupation_13_officeworker",
        "anyoccupation_14_healthcare", "anyoccupation_15_merchant",
        "anyoccupation_16_teacher", "anyoccupation_17_driver",
        "anyoccupation_18_sewerworker", "anyoccupation_19_guard",         
        "anyoccupation_20_unemployed", "anyoccupation_21_preworkingage",
        # Raw milk, blood, meat consumption,
        "pasteur_milk_py", "raw_milk_py",
        "raw_cow_blood_py","raw_goat_blood_py", "raw_sheep_blood_py",
        "raw_cow_py","raw_goat_py", "raw_sheep_py",
        # animal handling (milking, birthing, placenta, abortion,  waste disposal, room sharing)
        # "milkedanimals_pm", - this var not present Jan 2023
        "milkedcattle_py", "milkedcattle_pm", "milkedgoats_py", "milkedgoats_pm", "milkedsheep_py", "milkedsheep_pm",
        "birthing_py", "birthcattle_py", "birthcattle_pm", "birthgoats_py", "birthgoats_pm", "birthsheep_py", "birthsheep_pm",
        "handledplacenta_py", "placcattle_py", "placcattle_pm", "placgoats_py", "placgoats_pm", "placsheep_py", "placsheep_pm",
        "abortedproduct_py", "abortcattle_py", "abortcattle_pm", "abortgoats_py", "abortgoats_pm", "abortsheep_py",
        "animalwaste_py", "wastecattle_py", "wastecattle_pm", "wastegoats_py", "wastegoats_pm", "wastesheep_py", "wastesheep_pm",
        "sleepanimals_py", "sleepcattle_py", "sleepcattle_pm", "sleepgoats_py", "sleepgoats_pm", "sleepsheep_py", "sleepsheep_pm",
        "foeplac_py",
        # clinical symptoms (fever, lesion, rash, resp)
        "resp2wks","fever2wks", "diarrhoea2wks",
        "FeverTreatment",
        "rash_2weeks", "rash_py",
        "jointpain_2weeks", "jointpain_py",
        "cough_2weeks", "cough_py",
        "nightsweats_2weeks", "nightsweats_py",
        "backpain_2weeks", "backpain_py",
        # anthropometry
        "Final_weight", "Final_height", "WFH_zscore" 
        )]
# 

# merge H3 and P keeping all P data
HP<-merge(H3,P,by="HHID")
setdiff(P$HHID,H2$HHID) # all merges fine

# select fields for RVF serology outcomes

HS<-HS[,c("IndividualID","SN","rvfv")]
# 
PS<-merge(HP,HS,by="IndividualID")
# check the merge - 213 of 287 P records (and 296 HS records) merge
setdiff(HP$IndividualID,PS$IndividualID) # no serology for these individuals - check - possibly re aliquot n
table(Sam$species,Sam$SamType) # no reason to expect drop off in n samples available based on n collected
# n obs in HS with BZA prefix
summary(as.factor(substr(HS$IndividualID,1,3)))
# n = 222 BZA samples with results - 9 that do not merge given n results generated
setdiff(HS$IndividualID,PS$IndividualID) # all SEEDZ (makes sense) 
# from BacZoo - pilot individuals - no samples - n=9
# no mismatches - all generated serol data merge fine
##############################################################

# Keep human data separate for now

###################################################################################
###################################################################################
# summaries of outcome in rvf by species and ward

# humans
table(PS$WardID, PS$rvfv)
table(PS$Ward_Name, PS$rvfv)
table(PS$Village, PS$rvfv)
# seropos from 7 of 20 villages
table(PS$setting, PS$rvfv)
# n=13 of 15 seropos are from PA context
# raw prev in people
summary(as.factor(PS$rvfv))
binom.confint(15,15+198,method="exact")

# null model for clustering
null<-glmer(rvfv~ 1 + (1|HHID) + (1|WardID), 
              data=PS, family = binomial,
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=100000)))  
summary(null)
tab_model(null)

# animals
table(ARVF$setting, ARVF$rvfv_result)
# n=1 exposed animal in PU
table(ARVF$WardID, ARVF$rvfv_result)
table(ARVF$Ward_Name, ARVF$rvfv_result)

table(ARVF$Ward_Name, ARVF$rvfv_result, ARVF$species.x)

# overall by specise
table(ARVF$rvfv_result, ARVF$species.x, useNA="always")


###################################################################################
###################################################################################

# rm coord fields for households  (retain ward level)
library(dplyr)
ARVF<-select(ARVF,-c(CoordEW_Out, CoordNS_Out))
PS<-select(PS,-c(CoordEW_Out, CoordNS_Out))
# export csv files for Deng

# write.csv(ARVF,"BacZoo_Animals&RVFData_28Jan2023.csv")
# write.csv(PS,"BacZoo_Humans&RVFData_28Jan2023.csv")

###################################################################################
###################################################################################

# read in combined BZ & SZ data for livestock and people

BSH<-read.csv("/Users/mojohalliday/Library/CloudStorage/OneDrive-UniversityofGlasgow/General - Zoonoses Tanzania Data/Zoonoses Tanzania/SEEDZ/Final Data Versions/BZ and SEEDZ All serological results/RVF/rvf_human_data_2719.csv")
BSA<-read.csv("/Users/mojohalliday/Library/CloudStorage/OneDrive-UniversityofGlasgow/General - Zoonoses Tanzania Data/Zoonoses Tanzania/SEEDZ/Final Data Versions/BZ and SEEDZ All serological results/RVF/rvf_livestock_data_2719.csv")

# tables
summary(as.factor(BSA$species))
summary(as.factor(BSA$result))
table(BSA$species, BSA$result, BSA$study, useNA="always")
#

table(BSA$result[BSA$study=="baczoo"], 
      BSA$ward[BSA$study=="baczoo"], useNA="always")
