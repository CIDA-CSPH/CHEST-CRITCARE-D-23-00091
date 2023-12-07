################################################### -
## Title: Multiply impute missing data
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## PI: Ellen Burnham
## Date Created: 2023-11-01
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(mice, exclude="filter")
library(CIDAtools)

set.seed(2022)

# read analysis dataset
df_analysis <- read.csv(ProjectLocation("DataProcessed/merged_analysis_dataset.csv"))

##### Multiple imputations by chained equations (MICE) #####

# select relevant model variables
df_model <- df_analysis %>%
  select(record_id, year_enrolled,
         ###demographics
         age, sex_male, hispanic, race_nonwhite, body_mass_index,
         ###admission
         apache_2_admit, diabetes_y_n, cirrhosis_y_n, current_smoker, covid_cat,
         ###alcohol use
         peth, audit_c_patient, audit_c_nurse_validated, likely_current_drinking,
         ###continuous outcomes
         number_mv_days:hosp_free_days,
         ###dichotomous outcomes
         patient_discharged_to_home, patient_discharged_to_morge)

# initalize mice to obtain predictor matrix
init <- mice(df_model, defaultMethod=rep("", 4), maxit=0)

# modify predictor matrix to remove 'separated' data
pred_mat <- init$predictorMatrix
pred_mat["audit_c_patient", "patient_discharged_to_morge"] <- 0

# run mice to generate 5 imputed datasets
imp5 <- mice(
  df_model,
  defaultMethod=c("pmm", "logreg", "polyreg", ""),
  m=5, maxit=10,
  predictorMatrix=pred_mat,
  print=F)

# examine convergence of imputed peth, audit scores, and apache2 scores
plot(
  imp5,
  peth + audit_c_patient + audit_c_nurse_validated + apache_2_admit ~ .it | .ms,
  layout=c(2, 4))

# convert to long format for additional derivation
df_mice_long <- complete(imp5, include=T, action="long")

# stratify by pre/peri-pandemic
df_mice_long <- df_mice_long %>%
  mutate(pandemic_phase=factor(case_when(
    covid_cat=="Pre-pandemic" ~ "Pre-pandemic",
    TRUE ~ "Peri-pandemic"),
    c("Pre-pandemic", "Peri-pandemic")))

# categorize smoking history
df_mice_long <- df_mice_long %>%
  mutate(smoking_history=case_when(
    current_smoker==1 ~ "Current smoker",
    current_smoker==0 ~ "Never or former smoker"))

# bin PEth measures according to AFSHAR ET AL. (2022)
df_mice_long <- df_mice_long %>%
  mutate(peth_binned=factor(case_when(
    peth < 25 ~ "Low",
    peth >= 25 ~ "High"),
    c("Low", "High")))

# bin PEth measures according to ULWELLING ET AL. (2018)
df_mice_long <- df_mice_long %>%
  mutate(peth_binned_ulwelling=factor(case_when(
    peth < 20 ~ "Light/none",
    peth < 200 ~ "Significant",
    peth >= 200 ~ "Heavy"),
    c("Light/none", "Significant", "Heavy")))

# derive AUD based on RN-validated AUDIT-C
df_mice_long <- df_mice_long %>%
  mutate(any_aud_audit_c=case_when(
    audit_c_nurse_validated >= 4 & sex_male ~ 1,
    audit_c_nurse_validated >= 3 & !sex_male ~ 1,
    audit_c_nurse_validated < 5 ~ 0))

# compute log-transformed PEth
df_mice_long$peth_logged <- log(df_mice_long$peth + 10)

# convert back to multiply-imputed dataset object
imp5 <- as.mids(df_mice_long)

##### Write results to output file #####

# define output file name
out_name <- ProjectLocation("DataProcessed/multiply_imputed_dataset.rds")

# write to output file
saveRDS(imp5, file=out_name)
message("File written to: ", out_name)
