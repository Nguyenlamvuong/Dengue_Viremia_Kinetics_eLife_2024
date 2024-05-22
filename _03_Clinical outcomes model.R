#####################################################################################
# VIREMIA KINETICS
# Nguyen Lam Vuong
# 22 May 2024
# Landmark model for severe dengue & leakage outcome
#####################################################################################

# Load packages and functions
library(tidyverse)
library(rms)
library(geepack)

# Function to get predictions from gee model
get_pred_gee <- function(model, data) {
  datp <- data
  
  if (any(class(model) %in% "gee")) {
    tmp <- Epi::ci.exp(model, datp)
    datp$pred <- tmp[,1]/(1+tmp[,1])
    datp$low <- tmp[,2]/(1+tmp[,2])
    datp$upp <- tmp[,3]/(1+tmp[,3])
    
  } else {
    tmp <- predict(model, newdata = datp, type = "link", se.fit = T)
    critval <- qnorm(0.975) ## approx 95% CI
    upr <- tmp$fit + (critval * tmp$se.fit); upr2 <- model$family$linkinv(upr)
    lwr <- tmp$fit - (critval * tmp$se.fit); lwr2 <- model$family$linkinv(lwr)
    fit <- tmp$fit; fit2 <- model$family$linkinv(fit)
    datp$low <- lwr2 
    datp$upp <- upr2
    datp$pred <- fit2
  }
  
  return(datp)
}


# Function to get predictions from MI datasets
get_pred_mi <- function(model, newdata) {
  predm <- lapply(getfit(model), Epi::ci.lin, ctr.mat = newdata)
  Q <- U <- matrix(NA, nrow = nrow(newdata), ncol = length(predm))
  
  for (i in 1:length(predm)) {
    Q[, i] <- predm[[i]][, "Estimate"]
    U[, i] <- predm[[i]][, "StdErr"]^2
  }
  
  dfcom <- nrow(newdata) - 2
  
  ### pool predictions
  pred <- matrix(NA, nrow = nrow(Q), ncol = 2, dimnames = list(NULL, c("fit", "se.fit")))
  for(i in 1:nrow(Q)) {
    pi <- pool.scalar(Q[i, ], U[i, ], n = dfcom + 1)
    pred[i, 1] <- pi[["qbar"]]
    pred[i, 2] <- sqrt(pi[["t"]])
  }
  
  datp <- newdata
  pred_tmp <- exp(pred[,1])
  datp$pred <- pred_tmp/(1+pred_tmp)
  low_tmp <- exp(pred[,1] - qnorm(.975) * pred[,2])
  datp$low <- low_tmp/(1+low_tmp)
  upp_tmp <- exp(pred[,1] + qnorm(.975) * pred[,2])
  datp$upp <- upp_tmp/(1+upp_tmp)
  
  return(datp)
}


#####################################################################################
# Load data
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

## Data for landmark models
dat1 <- dat_vir %>%
  rename(Severe = Severe_dengue) %>%
  mutate(doi_severe = ifelse(Severe == 0, 10, doi_severe),
         doi_leakage = ifelse(is.na(doi_leakage) & Leakage == 1, 5, # Impute 5 for missing doi_leakage
                              ifelse(Leakage == 0, 10, doi_leakage))) # Max DOI of severe patients is 7

### Create landmark data, taking into account DOI of occurring severe dengue
dats1 <- dat1 %>% filter(doi_severe >= 1) %>% filter(DOI == 1) %>% mutate(LM = 1)
dats2 <- dat1 %>% filter(doi_severe >= 2) %>% filter(DOI == 2) %>% mutate(LM = 2)
dats3 <- dat1 %>% filter(doi_severe >= 3) %>% filter(DOI == 3) %>% mutate(LM = 3)
dats4 <- dat1 %>% filter(doi_severe >= 4) %>% filter(DOI == 4) %>% mutate(LM = 4)
dats5 <- dat1 %>% filter(doi_severe >= 5) %>% filter(DOI == 5) %>% mutate(LM = 5)
dats6 <- dat1 %>% filter(doi_severe >= 6) %>% filter(DOI == 6) %>% mutate(LM = 6)
dats7 <- dat1 %>% filter(doi_severe >= 7) %>% filter(DOI == 7) %>% mutate(LM = 7)

LMdats1 <- bind_rows(dats1, dats2, dats3, dats4, dats5, dats6, dats7) %>%
  arrange(Code, LM) %>%
  mutate(Code = factor(Code),
         cLM = LM - 3,
         clog10_vir = log10_vir - 6.38)

### Create landmark data for plasma leakage
datl1 <- dat1 %>% filter(doi_leakage >= 1) %>% filter(DOI == 1) %>% mutate(LM = 1)
datl2 <- dat1 %>% filter(doi_leakage >= 2) %>% filter(DOI == 2) %>% mutate(LM = 2)
datl3 <- dat1 %>% filter(doi_leakage >= 3) %>% filter(DOI == 3) %>% mutate(LM = 3)
datl4 <- dat1 %>% filter(doi_leakage >= 4) %>% filter(DOI == 4) %>% mutate(LM = 4)
datl5 <- dat1 %>% filter(doi_leakage >= 5) %>% filter(DOI == 5) %>% mutate(LM = 5)
datl6 <- dat1 %>% filter(doi_leakage >= 6) %>% filter(DOI == 6) %>% mutate(LM = 6)
datl7 <- dat1 %>% filter(doi_leakage >= 7) %>% filter(DOI == 7) %>% mutate(LM = 7)

### Combine to landmark dataset
LMdatl <- bind_rows(datl1, datl2, datl3, datl4, datl5, datl6, datl7) %>%
  arrange(Code, LM) %>%
  mutate(Code = factor(Code),
         cLM = LM - 3,
         clog10_vir = log10_vir - 6.38)


#####################################################################################
# Severe dengue-------------------------------------------------------------
## Complete-case analysis
### Fit model
LMdats1$Serology <- relevel(LMdats1$Serology, ref = "Probable secondary")

ms1 <- geeglm(Severe ~ rcs(clog10_vir, c(4,6.4,8.5)-6.38) + rcs(cLM,c(2,4,6)-3) + Serotype * Serology + rcs(cAge,c(7,13,28)-13) + Sex + vir01 + Study +
                clog10_vir : (rcs(cLM,c(2,4,6)-3) * Serotype * Serology + rcs(cAge,c(7,13,28)-13) + Sex + Study), 
              id = Code, data = LMdats1, family = binomial, corstr = "independence")

### Get predictions
nDF <- expand.grid(clog10_vir = c(log10(60), log10(300), log10(600), 3.01, 3:10) - 6.38,
                   cAge = c(5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30) - 13,
                   Sex = factor(c("Male","Female")),
                   cLM = c(1:7)-3,
                   Serotype = factor(c("DENV-1","DENV-2","DENV-3","DENV-4")),
                   Serology = factor(c("Probable primary","Probable secondary","Indeterminate")),
                   Study = factor(c("22DX", "DR", "MD"))) %>%
  mutate(log10_vir = clog10_vir + 6.38,
         clog10_vir2 = clog10_vir,
         vir01 = ifelse(log10_vir <= 3, 1, 0),
         Serology = factor(Serology, levels = c("Probable secondary","Probable primary","Indeterminate")),
         Age = cAge + 13,
         DOI = cLM + 3,
         Sex = factor(Sex, levels = c("Male", "Female")),
         Study = factor(Study, levels = c("22DX", "DR", "MD")),
         PCR = ifelse(Study == "22DX", "One-step", "Two-step"),
         PCR = factor(PCR, levels = c("Two-step", "One-step"))) %>%
  filter(!(clog10_vir<3-6.38 & PCR=="Two-step")) %>%
  filter(!(clog10_vir==log10(300)-6.38 & !Serotype %in% c("DENV-1", "DENV-3"))) %>%
  filter(!(clog10_vir==log10(60)-6.38 & !Serotype=="DENV-2")) %>%
  filter(!(clog10_vir==log10(600)-6.38 & !Serotype=="DENV-4")) %>%
  as.data.frame(.)

nDF1 <- nDF %>% filter(DOI==1) %>% filter(log10_vir>=5 & log10_vir<=10)
nDF2 <- nDF %>% filter(DOI==2) %>% filter((log10_vir>=4 & log10_vir<=10))
nDF3 <- nDF %>% filter(DOI==3) %>% filter((log10_vir>=3 & log10_vir<=9) | log10_vir<=3)
nDF4 <- nDF %>% filter(DOI==4) %>% filter((log10_vir>=3 & log10_vir<=8) | log10_vir<=3)
nDF5 <- nDF %>% filter(DOI==5) %>% filter((log10_vir>=3 & log10_vir<=7) | log10_vir<=3)
nDF6 <- nDF %>% filter(DOI==6) %>% filter((log10_vir>=3 & log10_vir<=6) | log10_vir<=3)
nDF7 <- nDF %>% filter(DOI==7) %>% filter((log10_vir>=3 & log10_vir<=5) | log10_vir<=3) %>% filter(Serotype!="DENV-4")

pms11 <- get_pred_gee(model = ms1, data = nDF1) %>% filter(!(log10_vir==3 & PCR=="One-step")) # to exclude the undetectable values setting as 1000 in 1-step PCR cohort, as this model used the LOD as LOD
pms12 <- get_pred_gee(model = ms1, data = nDF2) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pms13 <- get_pred_gee(model = ms1, data = nDF3) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pms14 <- get_pred_gee(model = ms1, data = nDF4) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pms15 <- get_pred_gee(model = ms1, data = nDF5) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pms16 <- get_pred_gee(model = ms1, data = nDF6) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pms17 <- get_pred_gee(model = ms1, data = nDF7) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pms1 <- bind_rows(pms11, pms12, pms13, pms14, pms15, pms16, pms17)


## Imputed data analysis
### Use imputed immune status from previous analysis
### 2 cases in this dataset are not in the dataset at enrollment (MD2628 and MD2695), but they all have immune status not-missing
library(mice)
load("imputation_sensitivity.RData")

all_dats <- complete(imp_enrol, action="long", include=T)

working_dats <- list()

### Use Serology from imputed data for landmark data
for(i in 0:max(all_dats$.imp)) {
  working_dats[[i+1]] <- all_dats %>%
    subset(.imp == i) %>%
    select(.imp, .id, Code, Serology) %>%
    left_join(select(LMdats1, -Serology), ., by = "Code") %>%
    mutate(.imp = i, .id = 1:nrow(.),
           Serology = ifelse(Code=="MD2628", "Primary",
                             ifelse(Code=="MD2695", "Secondary", as.character(Serology))),
           Serology = factor(Serology, levels = c("Secondary","Primary")),
           Code = factor(Code)) %>%
    relocate(.imp, .id, )
}

imputed_long <- as.mids(do.call(rbind, working_dats))

### Fit model
ms1i <- with(imputed_long,
             geeglm(Severe ~ rcs(clog10_vir, c(4,6.4,8.5)-6.38) + rcs(cLM,c(2,4,6)-3) + Serotype * Serology + rcs(cAge,c(7,13,28)-13) + Sex + vir01 + Study +
                      clog10_vir : (rcs(cLM,c(2,4,6)-3) * Serotype * Serology + rcs(cAge,c(7,13,28)-13) + Sex + Study), 
                    id = Code, family = binomial, corstr = "independence"))

### Get predictions
nDF <- expand.grid(clog10_vir = c(log10(60), log10(300), log10(600), 3.01, 3:10) - 6.38,
                   cAge = c(5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30) - 13,
                   Sex = factor(c("Male","Female")),
                   cLM = c(1:7)-3,
                   Serotype = factor(c("DENV-1","DENV-2","DENV-3","DENV-4")),
                   Serology = factor(c("Primary","Secondary")),
                   Study = factor(c("22DX", "DR", "MD"))) %>%
  mutate(log10_vir = clog10_vir + 6.38,
         clog10_vir2 = clog10_vir,
         vir01 = ifelse(log10_vir <= 3, 1, 0),
         Serology = factor(Serology, levels = c("Secondary","Primary")),
         Age = cAge + 13,
         DOI = cLM + 3,
         Sex = factor(Sex, levels = c("Male", "Female")),
         Study = factor(Study, levels = c("22DX", "DR", "MD")),
         PCR = ifelse(Study == "22DX", "One-step", "Two-step"),
         PCR = factor(PCR, levels = c("Two-step", "One-step"))) %>%
  filter(!(clog10_vir<3-6.38 & PCR=="Two-step")) %>%
  filter(!(clog10_vir==log10(300)-6.38 & !Serotype %in% c("DENV-1", "DENV-3"))) %>%
  filter(!(clog10_vir==log10(60)-6.38 & !Serotype=="DENV-2")) %>%
  filter(!(clog10_vir==log10(600)-6.38 & !Serotype=="DENV-4")) %>%
  as.data.frame(.)

nDF1 <- nDF %>% filter(DOI==1) %>% filter(log10_vir>=5 & log10_vir<=10)
nDF2 <- nDF %>% filter(DOI==2) %>% filter((log10_vir>=4 & log10_vir<=10))
nDF3 <- nDF %>% filter(DOI==3) %>% filter((log10_vir>=3 & log10_vir<=9) | log10_vir<=3)
nDF4 <- nDF %>% filter(DOI==4) %>% filter((log10_vir>=3 & log10_vir<=8) | log10_vir<=3)
nDF5 <- nDF %>% filter(DOI==5) %>% filter((log10_vir>=3 & log10_vir<=7) | log10_vir<=3)
nDF6 <- nDF %>% filter(DOI==6) %>% filter((log10_vir>=3 & log10_vir<=6) | log10_vir<=3)
nDF7 <- nDF %>% filter(DOI==7) %>% filter((log10_vir>=3 & log10_vir<=5) | log10_vir<=3) %>% filter(Serotype!="DENV-4")

pmsi11 <- get_pred_mi(model = ms1i, newdata = nDF1)
pmsi12 <- get_pred_mi(model = ms1i, newdata = nDF2)
pmsi13 <- get_pred_mi(model = ms1i, newdata = nDF3)
pmsi14 <- get_pred_mi(model = ms1i, newdata = nDF4)
pmsi15 <- get_pred_mi(model = ms1i, newdata = nDF5)
pmsi16 <- get_pred_mi(model = ms1i, newdata = nDF6)
pmsi17 <- get_pred_mi(model = ms1i, newdata = nDF7)

pmsi1 <- bind_rows(pmsi11, pmsi12, pmsi13, pmsi14, pmsi15, pmsi16, pmsi17) %>%
  mutate(Serology = factor(Serology, levels=c("Secondary","Primary"), labels=c("Probable secondary","Probable primary"))) %>%
  filter(!(log10_vir==3 & PCR=="One-step")) # to exclude the undetectable values setting as 1000 in 1-step PCR cohort, as this model used the LOD as LOD


# Plasma leakage-------------------------------------------------------------
## Complete-case analysis
### Fit model
LMdatl$Serology <- relevel(LMdatl$Serology, ref = "Probable secondary")

ml1 <- geeglm(Leakage ~ rcs(clog10_vir,c(4,6.4,8.5)-6.38) + vir01 + rcs(cLM,c(2,4,6)-3) + rcs(cAge,c(7,13,28)-13) + Sex + Serotype * Serology + Study +
                clog10_vir : (rcs(cLM,c(2,4,6)-3) * Serotype * Serology + rcs(cAge,c(7,13,28)-13) + Sex + Study) +
                vir01 : (cLM + cAge + Sex + Serotype + Serology + Study),
              id = Code, data = LMdatl, family = binomial, corstr = "independence")

### Get predictions
nDF <- expand.grid(clog10_vir = c(log10(60), log10(300), log10(600), 3.01, 3:10) - 6.38,
                   cAge = c(5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30) - 13,
                   Sex = factor(c("Male","Female")),
                   cLM = c(1:7)-3,
                   Serotype = factor(c("DENV-1","DENV-2","DENV-3","DENV-4")),
                   Serology = factor(c("Probable primary","Probable secondary","Indeterminate")),
                   Study = factor(c("22DX", "DR", "MD"))) %>%
  mutate(log10_vir = clog10_vir + 6.38,
         clog10_vir2 = clog10_vir,
         vir01 = ifelse(log10_vir <= 3, 1, 0),
         Serology = factor(Serology, levels = c("Probable secondary","Probable primary","Indeterminate")),
         Age = cAge + 13,
         DOI = cLM + 3,
         Sex = factor(Sex, levels = c("Male", "Female")),
         Study = factor(Study, levels = c("22DX", "DR", "MD")),
         PCR = ifelse(Study == "22DX", "One-step", "Two-step"),
         PCR = factor(PCR, levels = c("Two-step", "One-step"))) %>%
  filter(!(clog10_vir<3-6.38 & PCR=="Two-step")) %>%
  filter(!(clog10_vir==log10(300)-6.38 & !Serotype %in% c("DENV-1", "DENV-3"))) %>%
  filter(!(clog10_vir==log10(60)-6.38 & !Serotype=="DENV-2")) %>%
  filter(!(clog10_vir==log10(600)-6.38 & !Serotype=="DENV-4")) %>%
  as.data.frame(.)

nDF1 <- nDF %>% filter(DOI==1) %>% filter(log10_vir>=5 & log10_vir<=10)
nDF2 <- nDF %>% filter(DOI==2) %>% filter((log10_vir>=4 & log10_vir<=10))
nDF3 <- nDF %>% filter(DOI==3) %>% filter((log10_vir>=3 & log10_vir<=9) | log10_vir<=3)
nDF4 <- nDF %>% filter(DOI==4) %>% filter((log10_vir>=3 & log10_vir<=8) | log10_vir<=3)
nDF5 <- nDF %>% filter(DOI==5) %>% filter((log10_vir>=3 & log10_vir<=7) | log10_vir<=3)
nDF6 <- nDF %>% filter(DOI==6) %>% filter((log10_vir>=3 & log10_vir<=6) | log10_vir<=3)
nDF7 <- nDF %>% filter(DOI==7) %>% filter((log10_vir>=3 & log10_vir<=5) | log10_vir<=3) %>% filter(Serotype!="DENV-4")

pml11 <- get_pred_gee(model = ml1, data = nDF1) %>% filter(!(log10_vir==3 & PCR=="One-step")) # to exclude the undetectable values setting as 1000 in 1-step PCR cohort, as this model used the LOD as LOD
pml12 <- get_pred_gee(model = ml1, data = nDF2) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pml13 <- get_pred_gee(model = ml1, data = nDF3) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pml14 <- get_pred_gee(model = ml1, data = nDF4) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pml15 <- get_pred_gee(model = ml1, data = nDF5) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pml16 <- get_pred_gee(model = ml1, data = nDF6) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pml17 <- get_pred_gee(model = ml1, data = nDF7) %>% filter(!(log10_vir==3 & PCR=="One-step"))
pml1 <- bind_rows(pml11, pml12, pml13, pml14, pml15, pml16, pml17)


## Imputed data analysis
### Use imputed immune status from previous analysis
### 2 cases in this dataset are not in the dataset at enrollment (MD2628 and MD2695), but they all have immune status not-missing
library(mice)
load("imputation_sensitivity.RData")

all_dats <- complete(imp_enrol, action="long", include=T)

working_dats <- list()

### Use Serology from imputed data for landmark data
for(i in 0:max(all_dats$.imp)) {
  tmp <- all_dats %>%
    subset(.imp == i) %>%
    select(.imp, .id, Code, Serology, SILeak) %>%
    left_join(select(dat1, -Serology, -Leakage), ., by = "Code") %>%
    mutate(.imp = i, .id = 1:nrow(.),
           Serology = ifelse(Code=="MD2628", "Primary",
                             ifelse(Code=="MD2695", "Secondary", as.character(Serology))),
           Serology = factor(Serology, levels = c("Secondary","Primary")),
           Leakage = ifelse(Code %in% c("MD2628", "MD2695"), 0, SILeak),
           doi_leakage = ifelse(Leakage==0 & is.na(doi_leakage), 10,
                                ifelse(Leakage==1 & is.na(doi_leakage), 7, doi_leakage))) %>% 
    # impute doi_leakage by 7 to make all imputed dataset the same number of rows
    select(-SILeak) %>%
    relocate(.imp, .id, )
  
  ## Create landmark data for plasma leakage
  datl1 <- tmp %>% filter(doi_leakage >= 1) %>% filter(DOI == 1) %>% mutate(LM = 1)
  datl2 <- tmp %>% filter(doi_leakage >= 2) %>% filter(DOI == 2) %>% mutate(LM = 2)
  datl3 <- tmp %>% filter(doi_leakage >= 3) %>% filter(DOI == 3) %>% mutate(LM = 3)
  datl4 <- tmp %>% filter(doi_leakage >= 4) %>% filter(DOI == 4) %>% mutate(LM = 4)
  datl5 <- tmp %>% filter(doi_leakage >= 5) %>% filter(DOI == 5) %>% mutate(LM = 5)
  datl6 <- tmp %>% filter(doi_leakage >= 6) %>% filter(DOI == 6) %>% mutate(LM = 6)
  datl7 <- tmp %>% filter(doi_leakage >= 7) %>% filter(DOI == 7) %>% mutate(LM = 7)
  
  ## Combine to landmark dataset
  working_dats[[i+1]] <- bind_rows(datl1, datl2, datl3, datl4, datl5, datl6, datl7) %>%
    arrange(Code, LM) %>%
    mutate(Code = factor(Code),
           cLM = LM - 3,
           clog10_vir = log10_vir - 6.38)
}

### To solve the problem that working_dats[[1]] (the original dataset) had fewer number of rows than the imputed datasets
tmp190 <- as.data.frame(matrix(NA, nrow=190, ncol=27))
names(tmp190) <- names(working_dats[[1]])
tmp190$.imp <- 0
tmp190$.id <- c(11569:11758)
working_dats[[1]] <- bind_rows(working_dats[[1]], tmp190)

imputed_long <- as.mids(do.call(rbind, working_dats))

## Fit model
ml1i <- with(imputed_long,
             geeglm(Leakage ~ rcs(clog10_vir,c(4,6.4,8.5)-6.38) + vir01 + rcs(cLM,c(2,4,6)-3) + rcs(cAge,c(7,13,28)-13) + Sex + Serotype * Serology + Study +
                      clog10_vir : (rcs(cLM,c(2,4,6)-3) * Serotype * Serology + rcs(cAge,c(7,13,28)-13) + Sex + Study) +
                      vir01 : (cLM + cAge + Sex + Serotype + Serology + Study), 
                    id = Code, family = binomial, corstr = "independence")) # 1 min

## Get predictions
nDF <- expand.grid(clog10_vir = c(log10(60), log10(300), log10(600), 3.01, 3:10) - 6.38,
                   cAge = c(5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30) - 13,
                   Sex = factor(c("Male","Female")),
                   cLM = c(1:7)-3,
                   Serotype = factor(c("DENV-1","DENV-2","DENV-3","DENV-4")),
                   Serology = factor(c("Primary","Secondary")),
                   Study = factor(c("22DX", "DR", "MD"))) %>%
  mutate(log10_vir = clog10_vir + 6.38,
         clog10_vir2 = clog10_vir,
         vir01 = ifelse(log10_vir <= 3, 1, 0),
         Serology = factor(Serology, levels = c("Secondary","Primary")),
         Age = cAge + 13,
         DOI = cLM + 3,
         Sex = factor(Sex, levels = c("Male", "Female")),
         Study = factor(Study, levels = c("22DX", "DR", "MD")),
         PCR = ifelse(Study == "22DX", "One-step", "Two-step"),
         PCR = factor(PCR, levels = c("Two-step", "One-step"))) %>%
  filter(!(clog10_vir<3-6.38 & PCR=="Two-step")) %>%
  filter(!(clog10_vir==log10(300)-6.38 & !Serotype %in% c("DENV-1", "DENV-3"))) %>%
  filter(!(clog10_vir==log10(60)-6.38 & !Serotype=="DENV-2")) %>%
  filter(!(clog10_vir==log10(600)-6.38 & !Serotype=="DENV-4")) %>%
  as.data.frame(.)

nDF1 <- nDF %>% filter(DOI==1) %>% filter(log10_vir>=5 & log10_vir<=10)
nDF2 <- nDF %>% filter(DOI==2) %>% filter((log10_vir>=4 & log10_vir<=10))
nDF3 <- nDF %>% filter(DOI==3) %>% filter((log10_vir>=3 & log10_vir<=9) | log10_vir<=3)
nDF4 <- nDF %>% filter(DOI==4) %>% filter((log10_vir>=3 & log10_vir<=8) | log10_vir<=3)
nDF5 <- nDF %>% filter(DOI==5) %>% filter((log10_vir>=3 & log10_vir<=7) | log10_vir<=3)
nDF6 <- nDF %>% filter(DOI==6) %>% filter((log10_vir>=3 & log10_vir<=6) | log10_vir<=3)
nDF7 <- nDF %>% filter(DOI==7) %>% filter((log10_vir>=3 & log10_vir<=5) | log10_vir<=3) %>% filter(Serotype!="DENV-4")

pmli11 <- get_pred_mi(model = ml1i, newdata = nDF1)
pmli12 <- get_pred_mi(model = ml1i, newdata = nDF2)
pmli13 <- get_pred_mi(model = ml1i, newdata = nDF3)
pmli14 <- get_pred_mi(model = ml1i, newdata = nDF4)
pmli15 <- get_pred_mi(model = ml1i, newdata = nDF5)
pmli16 <- get_pred_mi(model = ml1i, newdata = nDF6)
pmli17 <- get_pred_mi(model = ml1i, newdata = nDF7)

pmli1 <- bind_rows(pmli11, pmli12, pmli13, pmli14, pmli15, pmli16, pmli17) %>%
  mutate(Serology = factor(Serology, levels=c("Secondary","Primary"), labels=c("Probable secondary","Probable primary"))) %>%
  filter(!(log10_vir==3 & PCR=="One-step")) # to exclude the undetectable values setting as 1000 in 1-step PCR cohort, as this model used the LOD as LOD


#####################################################################################
# Save for later use
save(pms1, pmsi1, pml1, pmli1, file = "Viremia and outcomes 240522.RData")
