#####################################################################################
# VIREMIA KINETICS - ANALYSIS FOR PLATELET COUNT
# Nguyen Lam Vuong
# 22 May 2024
#####################################################################################

# Load packages and functions
library(tidyverse)
library(lme4)
library(optimx) # to modify lmerControl
library(rms)


# Load data-------------------------------------------------------------------------------------------
## Full data
dat <- read_csv("Longitudinal viremia and platelet 22May2024.csv") %>%
  mutate(Study = factor(Study, levels = c("DR", "MD", "22DX")),
         Sex = factor(Sex, levels = c("Male", "Female")),
         Serotype = factor(Serotype, levels = c("DENV-1", "DENV-2", "DENV-3", "DENV-4", "Mixed", "Unknown")),
         Serology = factor(Serology, levels = c("Probable primary", "Probable secondary", "Indeterminate")),
         PLT_sq = PLT^(1/4),
         cAge = Age - 13,
         PCR = ifelse(Study == "22DX", 1, 0),
         PCR = factor(PCR, levels = c(0,1), labels = c("Two-step", "One-step")),
         vir_tmp = case_when(vir==0 & Serotype %in% c("DENV-1", "DENV-3") & PCR=="One-step" ~ 300,
                             vir==0 & Serotype=="DENV-2" & PCR=="One-step" ~ 60,
                             vir==0 & Serotype=="DENV-4" & PCR=="One-step" ~ 600,
                             vir==0 & PCR=="Two-step" ~ 1000,
                             vir>0 ~ vir),
         vir_tmp = ifelse(vir<1000 & PCR=="Two-step", 1000, vir_tmp), # 63 values in the 2-step PCR cohort are <1000 but detectable --> set them to undetectable
         vir01_original = ifelse(vir==0, 1, 0),
         vir01 = ifelse(vir<1000 & PCR=="Two-step", 1, vir01_original),
         log10_vir = log10(vir_tmp), 
         vir_10root = vir_tmp^.1) %>%
  select(Study, Code, Age, cAge, Sex, Enrol_DOI, Serotype, Serology, PCR, DOI, PLT, PLT_sq, 
         log10_vir, vir_10root, vir01_original, vir01, Severe_dengue, doi_severe, Leakage, doi_leakage)

## Create landmark data
dat0 <- dat %>%
  filter(DOI<=10) %>% # most were undetectable before day 10
  filter(!is.na(log10_vir)) %>%
  filter(!is.na(Serotype) & !Serotype %in% c("Mixed","Unknown")) %>% # exclude patients with negative PCR or mixed serotype
  mutate(Code = as.factor(Code),
         Serotype = factor(Serotype, levels=c("DENV-1","DENV-2","DENV-3","DENV-4")),
         Serology = relevel(Serology, ref="Probable secondary"))

dat_vir1 <- dat0 %>% filter(DOI==1) %>% select(Code, log10_vir, vir01)
dat_vir2 <- dat0 %>% filter(DOI==2) %>% select(Code, log10_vir, vir01)
dat_vir3 <- dat0 %>% filter(DOI==3) %>% select(Code, log10_vir, vir01)
dat_vir4 <- dat0 %>% filter(DOI==4) %>% select(Code, log10_vir, vir01)
dat_vir5 <- dat0 %>% filter(DOI==5) %>% select(Code, log10_vir, vir01)
dat_vir6 <- dat0 %>% filter(DOI==6) %>% select(Code, log10_vir, vir01)
dat_vir7 <- dat0 %>% filter(DOI==7) %>% select(Code, log10_vir, vir01)

dat_lm0 <- dat0 %>% select(Study, Code, PLT, PLT_sq, DOI, Age, cAge, Sex, Serotype, Serology, PCR)

dat_lm1 <- dat_lm0 %>% filter(DOI>=1) %>% filter(Code %in% dat_vir1$Code) %>% left_join(., dat_vir1, by="Code") %>% mutate(LM = 1)
dat_lm2 <- dat_lm0 %>% filter(DOI>=2) %>% filter(Code %in% dat_vir2$Code) %>% left_join(., dat_vir2, by="Code") %>% mutate(LM = 2)
dat_lm3 <- dat_lm0 %>% filter(DOI>=3) %>% filter(Code %in% dat_vir3$Code) %>% left_join(., dat_vir3, by="Code") %>% mutate(LM = 3)
dat_lm4 <- dat_lm0 %>% filter(DOI>=4) %>% filter(Code %in% dat_vir4$Code) %>% left_join(., dat_vir4, by="Code") %>% mutate(LM = 4)
dat_lm5 <- dat_lm0 %>% filter(DOI>=5) %>% filter(Code %in% dat_vir5$Code) %>% left_join(., dat_vir5, by="Code") %>% mutate(LM = 5)
dat_lm6 <- dat_lm0 %>% filter(DOI>=6) %>% filter(Code %in% dat_vir6$Code) %>% left_join(., dat_vir6, by="Code") %>% mutate(LM = 6)
dat_lm7 <- dat_lm0 %>% filter(DOI>=7) %>% filter(Code %in% dat_vir7$Code) %>% left_join(., dat_vir7, by="Code") %>% mutate(LM = 7)

LMdat <- bind_rows(dat_lm1, dat_lm2, dat_lm3, dat_lm4, dat_lm5, dat_lm6, dat_lm7) %>% 
  arrange(Code, LM, DOI) %>% 
  filter(!is.na(PLT)) %>%
  mutate(LM_factor = as.factor(LM))


#####################################################################################
# Landmark supermodel
## Include the interaction between viremia and Serotype, Serology, and PCR
## Include the interaction between DOI and LM
## Random effect with splines of DOI
mplm8 <- lmer(PLT_sq ~ (rcs(log10_vir, c(4,6,8)) + vir01) * Serotype * (Serology + PCR) + rcs(cAge, c(-6,0,15)) + Sex + rcs(DOI, c(2,5,8)) + rcs(LM, c(2,4,6)) +
                rcs(DOI, c(2,5,8)) : ((rcs(log10_vir, c(4,6,8)) + vir01) * Serotype * PCR + rcs(cAge, c(-6,0,15)) + Sex + Serology + rcs(LM, c(2,4,6))) + 
                rcs(LM, c(2,4,6)) : ((rcs(log10_vir, c(4,6,8)) + vir01) * Serotype * PCR) +
                (rcs(DOI, c(2,5,8))|Code:LM), 
              REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')), data=LMdat) # 9 mins
# fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients


#------------------------------------
## Save for later use
save(mplm8, file = "Platelet landmark model 240522.RData")
