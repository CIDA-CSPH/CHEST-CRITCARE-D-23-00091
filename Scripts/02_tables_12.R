################################################### -
## Title: Create manuscript tables 1 and 2
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## PI: Ellen Burnham
## Date Created: 2023-11-01
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(gtsummary)
library(CIDAtools)

# read analysis dataset
df_analysis <- read.csv(ProjectLocation("DataProcessed/merged_analysis_dataset.csv"))

##### Create Table 1 #####

# Table 1: Study demographics and in-patient outcomes, stratified by pandemic era
tab_demog_outcomes_strat <- df_analysis %>%
  select(pandemic_phase, age, sex_male, hispanic, race,
         race_nonwhite, smoking_history, apache_2_admit,
         diabetes_y_n, cirrhosis_y_n, body_mass_index,
         number_mv_days, total_no_icu_days, total_number_of_days_in_ho,
         patient_discharged_to_home, patient_discharged_to_morge) %>%
  mutate(pandemic_phase=factor(pandemic_phase, c("Pre-pandemic", "Peri-pandemic"))) %>%
  tbl_summary(
    by=pandemic_phase,
    label=list(
      age ~ "Age, mean (SD)",
      sex_male ~ "Male sex",
      hispanic ~ "Hispanic ethnicity",
      race ~ "Race (self-reported)",
      race_nonwhite ~ "Non-white race",
      smoking_history ~ "Smoking history",
      apache_2_admit ~ "APACHE II score at admission, mean (SD)",
      diabetes_y_n ~ "History of Diabetes",
      cirrhosis_y_n ~ "History of Cirrhosis",
      body_mass_index ~ "Body Mass Index [kg/m2], mean (SD)",
      number_mv_days ~ "Days of mechanical ventilation, median (IQR)",
      total_no_icu_days ~ "Days in ICU, median (IQR)",
      total_number_of_days_in_ho ~ "Days in hospital, median (IQR)",
      patient_discharged_to_home ~ "Discharged to home",
      patient_discharged_to_morge ~ "In-hospital mortality"),
    statistic=list(
      age ~ "{mean} ({sd})",
      apache_2_admit ~ "{mean} ({sd})",
      body_mass_index ~ "{mean} ({sd})"),
    digits=list(
      age ~ c(1, 1),
      apache_2_admit ~ c(1, 1),
      body_mass_index ~ c(1, 1)),
    missing="no") %>%
  add_p(
    test=list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "fisher.test")) %>%
  bold_p() %>%
  modify_header(label="**Clinical Characteristics**") %>%
  modify_footnote(everything() ~ NA) %>%
  modify_caption(
    "Table 1: Study demographics and in-patient outcomes, stratified by pandemic era")

##### Create Table 2 #####

# Table 2: Alcohol use in the cohort
tab_aud_strat <- df_analysis %>%
  select(pandemic_phase, audit_c_nurse_validated, peth, any_aud_audit_c,
         peth_binned, likely_current_drinking, peth_binned_ulwelling) %>%
  mutate(pandemic_phase=factor(pandemic_phase, c("Pre-pandemic", "Peri-pandemic"))) %>%
  mutate(peth_binned_ulwelling=factor(
    peth_binned_ulwelling,
    c("Light/none", "Significant", "Heavy"),
    c("Light or no alcohol consumption (PEth <20 ng/mL)",
      "Significant alcohol consumption (PEth 20-199 ng/mL)",
      "Heavy alcohol consumption (PEth ≥ 200)"))) %>%
  mutate(audit_c_zero=case_when(
    audit_c_nurse_validated==0 ~ 1,
    audit_c_nurse_validated!=0 ~ 0)) %>%
  tbl_summary(
    by=pandemic_phase,
    label=list(
      audit_c_zero ~ "AUDIT-C = Zero",
      audit_c_nurse_validated ~ "AUDIT-C score, median (IQR)",
      peth ~ "PEth [ng/mL], median (IQR)",
      any_aud_audit_c ~ "Likely alcohol misuse, AUDIT-C definition",
      peth_binned ~ "Alcohol Misuse, PEth cut-point from AFSHAR ET AL.",
      likely_current_drinking ~ "Current alcohol drinking, composite definition",
      peth_binned_ulwelling ~ "Alcohol consumption, PEth definition"),
    missing="no") %>%
  add_p(
    test=list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "fisher.test")) %>%
  bold_p() %>%
  modify_header(label="**Patient Alcohol Use Metrics and Characteristics**") %>%
  modify_footnote(everything() ~ NA) %>%
  modify_caption("Table 2: Alcohol use in the cohort.")

# summarize PEth levels for those with AUDIT-C==0
df_zero_audit <- subset(df_analysis,
                        !is.na(audit_c_nurse_validated) & audit_c_nurse_validated==0)
summary(df_zero_audit$peth)
median_peth_zero_audit <- median(df_zero_audit$peth, na.rm=T)
median_peth_zero_audit_pre <- median(
  df_zero_audit$peth[df_zero_audit$prepandemic==1], na.rm=T)
median_peth_zero_audit_peri <- median(
  df_zero_audit$peth[df_zero_audit$peripandemic==1], na.rm=T)


# availability of AUDIT-C
audit_c_avail <- (!is.na(df_analysis$audit_c_patient) | 
                    !is.na(df_analysis$audit_c_nurse_validated))

audit_c_avail_pre <- mean(audit_c_avail[df_analysis$prepandemic == 1])
audit_c_avail_peri <- mean(audit_c_avail[df_analysis$prepandemic == 0])

# availability of PEth by cohort
peth_avail_pre <- 1 - mean(is.na(df_analysis$peth) & df_analysis$prepandemic)
peth_avail_peri <- 1 - mean(is.na(df_analysis$peth) & df_analysis$peripandemic)
