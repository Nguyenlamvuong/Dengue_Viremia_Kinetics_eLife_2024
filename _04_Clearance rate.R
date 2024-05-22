#####################################################################################
# VIREMIA KINETICS
# Nguyen Lam Vuong
# 22 May 2024
# Viral clearance rate (slope) and clinical outcomes
#####################################################################################

# Load packages and functions
library(tidyverse)
library(rms)
library(gtsummary)
library(MCMCglmm)


#####################################################################################
# Full data
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

# Data for viremia kinetics
dat_vir <- dat %>%
  filter(DOI<=10) %>% # most were undetectable before day 10
  filter(!is.na(log10_vir)) %>%
  filter(!is.na(Serotype) & !Serotype %in% c("Mixed","Unknown")) %>% # exclude patients with negative PCR
  mutate(Serotype = factor(Serotype, levels=c("DENV-1","DENV-2","DENV-3","DENV-4")))

dat_vir12 <- dat_vir %>%
  filter(!(Serotype == "DENV-4" & PCR == "Two-step")) %>%
  mutate(log10_vir_fake = ifelse(vir01==1, log10_vir-1, log10_vir),
         ymax = log10_vir,
         ymin = ifelse(vir01==1, -Inf, log10_vir)) %>%
  as.data.frame(.)


#####################################################################################
# Derive slope of viremia decline using the random effect of the model of viremia kinetics
## Fit models
### Treat illness day (DOI) as a linear effect rather than splines
m23 <- MCMCglmm(cbind(ymin, ymax) ~ DOI + rcs(cAge, c(-6,0,15)) + Sex + Serotype + Serology + PCR + 
                  Serotype : Serology + 
                  DOI : (rcs(cAge, c(-6,0,15)) + Sex + Serotype + Serology + PCR),
                random = ~ idh(1 + DOI) : Code,
                family = "cengaussian", nitt=20000, pr=T, data = dat_vir12) # 5 mins in Mac

## Extract random slope from model with linear trend of DOI in the random effect
library(broom.mixed)
random_data <- tidy(m23, effects="ran_vals") %>% as.data.frame()
random_slope <- random_data %>% filter(term=="DOI") %>% select(level, estimate)
random_intercept <- random_data %>% filter(term=="(Intercept)") %>% select(level, estimate)
random_indi <- full_join(random_slope, random_intercept, by = "level")
names(random_indi) <- c("Code", "Random_slope", "Random_intercept")

## Predict using the same data
pred23_original <- predict(m23) %>% as.data.frame(.)

## Data for prediction to get slope of individuals
nDF1 <- dat_vir12 %>%
  group_by(Code) %>%
  slice(1) %>%
  ungroup() %>%
  select(Code, cAge, Sex, Serotype, Serology, PCR)

nDF2 <- expand.grid(Code = nDF1$Code, DOI = 2:5)

nDF <- left_join(nDF2, nDF1, by = "Code") %>%
  arrange(Code, DOI) %>%
  mutate(ymin = 0,
         ymax = 10)

## DOI with linear trend in the random effect
pred23 <- predict(m23, newdata = nDF) %>% as.data.frame(.) %>% rename(fit = V1)
pm23 <- bind_cols(nDF, pred23)

ran_slope23 <- pm23 %>%
  left_join(., random_indi, by = "Code") %>%
  mutate(newfit = fit + Random_slope * DOI + Random_intercept) %>%
  select(Code, DOI, newfit) %>%
  pivot_wider(id_cols = Code, names_from = DOI, names_prefix = "doi", names_sep = "", values_from = newfit) %>%
  mutate(Ran_slope_linear = doi2 - doi3) %>%
  select(Code, Ran_slope_linear)

## Combine data
dat_slope <- dat_vir12 %>%
  group_by(Code) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(., ran_slope23, by = "Code") %>%
  mutate(Slope = Ran_slope_linear) %>%
  select(Code, Study, Age, Sex, Enrol_DOI, Serotype, Serology, PCR, Slope, Severe_dengue, Leakage)


#####################################################################################
# Save data with viremia decline rate------------------------------------------------
save(dat_slope, file = "Data with viremia decline rate 240522.RData")
