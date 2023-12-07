################################################### -
## Title: Create tables for supplementary materials
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## PI: Ellen Burnham
## Date Created: 2023-11-01
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(gtsummary)
library(mice, exclude="filter")
library(CIDAtools)

# load multiply-imputed data
imp5 <- readRDS(ProjectLocation("DataProcessed/multiply_imputed_dataset.rds"))

# isolate original (unimputed) data 
df_unimputed <- complete(imp5, 0)

##### Helper functions (copy from 03_table_3.R) #####

# define helper function to pool ORM results for a given effect
my_orm_pooler <- function(fit, eff_name, k=1) {
  estimates <- unlist(lapply(fit$analyses, function(x) x$coefficients[eff_name]))
  vars <- unlist(lapply(fit$analyses, function(x) x$var[eff_name, eff_name]))
  n <- unname(fit$analyses[[1]]$stats["Obs"])
  res <- pool.scalar(Q=estimates, U=vars, n=n, k=k)
  t_stat <- res$qbar / sqrt(res$t)
  p_val <- pt(abs(t_stat), df=res$df, lower.tail=F) * 2
  t_crit <- rep(qt(p=0.025, df=res$df, lower.tail=F), 2)
  t_crit[1] <- -t_crit[1] ###for lower bound
  ci_95 <- res$qbar + sqrt(res$t) * t_crit
  return(list(coef=res$qbar, std.err=sqrt(res$t), df=res$df,
              ci_lwr=ci_95[1], ci_upr=ci_95[2],
              t.stat=t_stat, p.value=p_val))
}

# define helper function to aggregate list of effects for ORM results
my_pool_summary <- function(fit, eff_names, eff_labels=NULL,  k=1){
  if(length(eff_labels) < length(eff_names)) eff_labels <- eff_names
  df_out <- data.frame(
    effect=eff_labels[1], my_orm_pooler(fit, eff_name=eff_names[1], k=k))
  if(length(eff_labels) > 1){
    for(i in 2:length(eff_names)){
      df_out <- rbind(
        df_out,
        data.frame(
          effect=eff_labels[i], my_orm_pooler(fit, eff_name=eff_names[i], k=k)))
    }    
  }
  return(df_out)
}

##### Supplemental Table 1: Relationship between PEth and covariates #####

# perform pooled analysis with PEth as an outcome
fit_peth_covar <- with(
  imp5,
  rms::orm(peth ~ age + sex_male + race_nonwhite + hispanic + apache_2_admit +
             diabetes_y_n + cirrhosis_y_n + current_smoker + covid_cat))
tbl_peth_covar <- my_pool_summary(
  fit_peth_covar,
  c("age", "sex_male", "race_nonwhite", "hispanic", "apache_2_admit",
    "diabetes_y_n", "cirrhosis_y_n", "current_smoker",
    "covid_cat=Positive", "covid_cat=Pre-pandemic"),
  c("Age", "Sex: Male", "Race: Non-white", "Hispanic", "APACHE II",
    "Diabetes +", "Cirrhosis +", "Current smoker",
    "COVID +", "Pre-pandemic")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

##### Supplemental Table 2: effects of PEth (trichotomized) on selected outcomes #####

# perform pooled analysis with days on MV as an outcome
fit_mv_days_unadj <- with(imp5, rms::orm(number_mv_days ~ peth_binned_ulwelling))
tbl_mv_days_unadj <- my_pool_summary(
  fit_mv_days_unadj,
  c("peth_binned_ulwelling=Significant", "peth_binned_ulwelling=Heavy")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with ICU days as an outcome
fit_icu_days_unadj <- with(imp5, rms::orm(total_no_icu_days ~ peth_binned_ulwelling))
tbl_icu_days_unadj <- my_pool_summary(
  fit_icu_days_unadj,
  c("peth_binned_ulwelling=Significant", "peth_binned_ulwelling=Heavy")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with hospital days as an outcome
fit_hosp_days_unadj <- with(imp5, rms::orm(total_number_of_days_in_ho ~ peth_binned_ulwelling))
tbl_hosp_days_unadj <- my_pool_summary(
  fit_hosp_days_unadj,
  c("peth_binned_ulwelling=Significant", "peth_binned_ulwelling=Heavy")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with in-hospital mortality as an outcome
fit_mortality_unadj <- with(imp5, rms::orm(patient_discharged_to_morge ~ peth_binned_ulwelling))
tbl_mortality_unadj <- my_pool_summary(
  fit_mortality_unadj,
  c("peth_binned_ulwelling=Significant", "peth_binned_ulwelling=Heavy")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with discharge as an outcome
fit_discharge_unadj <- with(imp5, rms::orm(patient_discharged_to_home ~ peth_binned_ulwelling))
tbl_discharge_unadj <- my_pool_summary(
  fit_discharge_unadj,
  c("peth_binned_ulwelling=Significant", "peth_binned_ulwelling=Heavy")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# for the following outcomes all of the predictors are identical
eff_names <-c("age", "sex_male", "race_nonwhite", "hispanic", "apache_2_admit",
              "diabetes_y_n", "cirrhosis_y_n", "current_smoker",
              "covid_cat=Positive", "covid_cat=Pre-pandemic",
              "peth_binned_ulwelling=Significant", "peth_binned_ulwelling=Heavy")
eff_labels <- c("Age", "Sex: Male", "Race: Non-white", "Hispanic", "APACHE II",
                "Diabetes +", "Cirrhosis +", "Current smoker",
                "COVID +", "Pre-pandemic", "PEth level: Significant",
                "PEth level: Heavy")

# perform pooled analysis with days on MV as an outcome
fit_mv_days_adj <- with(
  imp5,
  rms::orm(number_mv_days ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned_ulwelling))
tbl_mv_days_adj <- my_pool_summary(fit_mv_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with ICU days as an outcome
fit_icu_days_adj <- with(
  imp5,
  rms::orm(total_no_icu_days ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned_ulwelling))
tbl_icu_days_adj <- my_pool_summary(fit_icu_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with hospital days as an outcome
fit_hosp_days_adj <- with(
  imp5,
  rms::orm(total_number_of_days_in_ho ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned_ulwelling))
tbl_hosp_days_adj <- my_pool_summary(fit_hosp_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with in-hospital mortality as an outcome
fit_mortality_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_morge ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned_ulwelling))
tbl_mortality_adj <- my_pool_summary(fit_mortality_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with discharge as an outcome
fit_discharge_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_home ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned_ulwelling))
tbl_discharge_adj <- my_pool_summary(fit_discharge_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# begin tabulating unadjusted effects of PEth
tbl_mv_days_unadj$Outcome <- "Mechanical ventilation days"
tbl_icu_days_unadj$Outcome <- "Intensive care unit days"
tbl_hosp_days_unadj$Outcome <- "Hospital days"
tbl_mortality_unadj$Outcome <- "Mortality"
tbl_discharge_unadj$Outcome <- "Discharge to home"
tbl_s2_unadj <- bind_rows(
  tbl_mv_days_unadj,
  tbl_icu_days_unadj,
  tbl_hosp_days_unadj,
  tbl_mortality_unadj,
  tbl_discharge_unadj) %>%
  mutate(effect=case_when(
    effect=="peth_binned_ulwelling=Significant" ~ "PEth level: Significant",
    effect=="peth_binned_ulwelling=Heavy" ~ "PEth level: Heavy")) %>%
  select(Outcome, effect, `odds ratio`, contains("95%"), `p-value`)

# begin tabulating adjusted effects of PEth
tbl_mv_days_adj$Outcome <- "Mechanical ventilation days"
tbl_icu_days_adj$Outcome <- "Intensive care unit days"
tbl_hosp_days_adj$Outcome <- "Hospital days"
tbl_mortality_adj$Outcome <- "Mortality"
tbl_discharge_adj$Outcome <- "Discharge to home"
tbl_s2_adj <- bind_rows(
  filter(tbl_mv_days_adj, grepl("PEth level", effect)),
  filter(tbl_icu_days_adj, grepl("PEth level", effect)),
  filter(tbl_hosp_days_adj, grepl("PEth level", effect)),
  filter(tbl_mortality_adj, grepl("PEth level", effect)),
  filter(tbl_discharge_adj, grepl("PEth level", effect))) %>%
  select(Outcome, effect, `odds ratio`, contains("95%"), `p-value`)

##### Supplemental Table 3: effects of AUDIT-C on selected outcomes #####

# perform pooled analysis with days on MV as an outcome
fit_mv_days_unadj <- with(imp5, rms::orm(number_mv_days ~ any_aud_audit_c))
tbl_mv_days_unadj <- my_pool_summary(fit_mv_days_unadj, c("any_aud_audit_c")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with ICU days as an outcome
fit_icu_days_unadj <- with(imp5, rms::orm(total_no_icu_days ~ any_aud_audit_c))
tbl_icu_days_unadj <- my_pool_summary(fit_icu_days_unadj, c("any_aud_audit_c")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with hospital days as an outcome
fit_hosp_days_unadj <- with(imp5, rms::orm(total_number_of_days_in_ho ~ any_aud_audit_c))
tbl_hosp_days_unadj <- my_pool_summary(fit_hosp_days_unadj, c("any_aud_audit_c")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with in-hospital mortality as an outcome
fit_mortality_unadj <- with(imp5, rms::orm(patient_discharged_to_morge ~ any_aud_audit_c))
tbl_mortality_unadj <- my_pool_summary(fit_mortality_unadj, c("any_aud_audit_c")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with discharge as an outcome
fit_discharge_unadj <- with(imp5, rms::orm(patient_discharged_to_home ~ any_aud_audit_c))
tbl_discharge_unadj <- my_pool_summary(fit_discharge_unadj, c("any_aud_audit_c")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# for the following outcomes all of the predictors are identical
eff_names <-c("age", "sex_male", "race_nonwhite", "hispanic", "apache_2_admit",
              "diabetes_y_n", "cirrhosis_y_n", "current_smoker",
              "covid_cat=Positive", "covid_cat=Pre-pandemic",
              "any_aud_audit_c")
eff_labels <- c("Age", "Sex: Male", "Race: Non-white", "Hispanic", "APACHE II",
                "Diabetes +", "Cirrhosis +", "Current smoker",
                "COVID +", "Pre-pandemic", "Likely Alc. Misuse, AUDIT-C")

# perform pooled analysis with days on MV as an outcome
fit_mv_days_adj <- with(
  imp5,
  rms::orm(number_mv_days ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + any_aud_audit_c))
tbl_mv_days_adj <- my_pool_summary(fit_mv_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with ICU days as an outcome
fit_icu_days_adj <- with(
  imp5,
  rms::orm(total_no_icu_days ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + any_aud_audit_c))
tbl_icu_days_adj <- my_pool_summary(fit_icu_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with hospital days as an outcome
fit_hosp_days_adj <- with(
  imp5,
  rms::orm(total_number_of_days_in_ho ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + any_aud_audit_c))
tbl_hosp_days_adj <- my_pool_summary(fit_hosp_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with in-hospital mortality as an outcome
fit_mortality_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_morge ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + any_aud_audit_c))
tbl_mortality_adj <- my_pool_summary(fit_mortality_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with discharge as an outcome
fit_discharge_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_home ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + any_aud_audit_c))
tbl_discharge_adj <- my_pool_summary(fit_discharge_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# begin tabulating unadjusted effects of AUDIT-C
tbl_mv_days_unadj$Outcome <- "Mechanical ventilation days"
tbl_icu_days_unadj$Outcome <- "Intensive care unit days"
tbl_hosp_days_unadj$Outcome <- "Hospital days"
tbl_mortality_unadj$Outcome <- "Mortality"
tbl_discharge_unadj$Outcome <- "Discharge to home"
tbl_s3_unadj <- bind_rows(
  tbl_mv_days_unadj,
  tbl_icu_days_unadj,
  tbl_hosp_days_unadj,
  tbl_mortality_unadj,
  tbl_discharge_unadj) %>%
  select(Outcome, `odds ratio`, contains("95%"), `p-value`)

# begin tabulating adjusted effects of AUDIT-C
tbl_mv_days_adj$Outcome <- "Mechanical ventilation days"
tbl_icu_days_adj$Outcome <- "Intensive care unit days"
tbl_hosp_days_adj$Outcome <- "Hospital days"
tbl_mortality_adj$Outcome <- "Mortality"
tbl_discharge_adj$Outcome <- "Discharge to home"
tbl_s3_adj <- bind_rows(
  filter(tbl_mv_days_adj, effect=="Likely Alc. Misuse, AUDIT-C"),
  filter(tbl_icu_days_adj, effect=="Likely Alc. Misuse, AUDIT-C"),
  filter(tbl_hosp_days_adj, effect=="Likely Alc. Misuse, AUDIT-C"),
  filter(tbl_mortality_adj, effect=="Likely Alc. Misuse, AUDIT-C"),
  filter(tbl_discharge_adj, effect=="Likely Alc. Misuse, AUDIT-C")) %>%
  select(Outcome, `odds ratio`, contains("95%"), `p-value`)


##### Supplemental Table 4: Relationship between PEth and MV Days #####

# for the following outcomes all of the predictors are identical
eff_names <-c("age", "sex_male", "race_nonwhite", "hispanic", "apache_2_admit",
              "diabetes_y_n", "cirrhosis_y_n", "current_smoker",
              "covid_cat=Positive", "covid_cat=Pre-pandemic",
              "peth_binned=High")
eff_labels <- c("Age", "Sex: Male", "Race: Non-white", "Hispanic", "APACHE II",
                "Diabetes +", "Cirrhosis +", "Current smoker",
                "COVID +", "Pre-pandemic", "PEth level: High")

# perform pooled analysis with days on MV as an outcome
fit_mv_days_adj <- with(
  imp5,
  rms::orm(number_mv_days ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_s4 <- my_pool_summary(fit_mv_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

##### Supplemental Table 5: Relationship between PEth and Days in ICU #####

# perform pooled analysis with ICU days as an outcome
fit_icu_days_adj <- with(
  imp5,
  rms::orm(total_no_icu_days ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_s5 <- my_pool_summary(fit_icu_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

##### Supplemental Table 6: Relationship between PEth and Mortality  #####

# perform pooled analysis with in-hospital mortality as an outcome
fit_mortality_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_morge ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_s6 <- my_pool_summary(fit_mortality_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)


##### Supplemental Table 7: Relationship between PEth and Home Discharge #####

# perform pooled analysis with discharge as an outcome
fit_discharge_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_home ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_s7 <- my_pool_summary(fit_discharge_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)
