#####################################################################################
# VIREMIA KINETICS
# Nguyen Lam Vuong
# 22 May 2024
# Viremia models with censored data
#####################################################################################

# Load packages and functions
library(tidyverse)
library(rms)
library(MCMCglmm)


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

## Data for viremia kinetics
dat_vir <- dat %>%
  filter(DOI<=10) %>% # most were undetectable before day 10
  filter(!is.na(log10_vir)) %>%
  filter(!is.na(Serotype) & !Serotype %in% c("Mixed","Unknown")) %>% # exclude patients with negative PCR
  mutate(Serotype = factor(Serotype, levels=c("DENV-1","DENV-2","DENV-3","DENV-4")))


# Data for log-10 viremia model------------------------------------------------------------------------
dat_vir12 <- dat_vir %>%
  filter(!(Serotype == "DENV-4" & PCR == "Two-step")) %>%
  mutate(log10_vir_fake = ifelse(vir01==1, log10_vir-1, log10_vir),
         ymax = log10_vir,
         ymin = ifelse(vir01==1, -Inf, log10_vir)) %>%
  as.data.frame(.)


# Data for 10-root viremia model------------------------------------------------------------------------
dat_vir12b <- dat_vir %>%
  filter(!(Serotype == "DENV-4" & PCR == "Two-step")) %>%
  mutate(vir_10root_fake = ifelse(vir01==1, vir_10root-1, vir_10root),
         ymax = vir_10root,
         ymin = ifelse(vir01==1, 0, vir_10root)) %>%
  as.data.frame(.)


#####################################################################################
# Model
## Log-10 viremia------------------------------------
m18 <- MCMCglmm(cbind(ymin, ymax) ~ rcs(DOI, c(1,2,4,6)) + rcs(cAge, c(-6,0,15)) + Sex + Serotype + Serology + PCR + 
                  Serotype : Serology + 
                  rcs(DOI, c(1,2,4,6)) : (rcs(cAge, c(-6,0,15)) + Sex + Serotype + Serology + PCR),
                random = ~ idh(1 + rcs(DOI, c(1,2,4,6))) : Code,
                family = "cengaussian", nitt=20000, data = dat_vir12) # 8 mins (OUCRU PC)

## 10-root viremia------------------------------------
m18b <- MCMCglmm(cbind(ymin, ymax) ~ rcs(DOI, c(1,2,4,6)) + rcs(cAge, c(-6,0,15)) + Sex + Serotype + Serology + PCR + 
                   Serotype : Serology + 
                   rcs(DOI, c(1,2,4,6)) : (rcs(cAge, c(-6,0,15)) + Sex + Serotype + Serology + PCR),
                 random = ~ idh(1 + rcs(DOI, c(1,2,4,6))) : Code,
                 family = "cengaussian", nitt=20000, data = dat_vir12b) # 8 mins (OUCRU PC)


#####################################################################################
# Save models--------------------------------------------------------
save(m18, m18b, file = "Viremia model censored 240522.RData")

