################################################### -
## Title: Create manuscript table 3
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

##### Helper functions #####

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


##### Fit unadjusted ORMs to multiply imputed datasets #####

# perform pooled analysis with days on MV as an outcome
fit_mv_days_unadj <- with(imp5, rms::orm(number_mv_days ~ peth_binned))
tbl_mv_days_unadj <- my_pool_summary(fit_mv_days_unadj, c("peth_binned=High")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with ICU days as an outcome
fit_icu_days_unadj <- with(imp5, rms::orm(total_no_icu_days ~ peth_binned))
tbl_icu_days_unadj <- my_pool_summary(fit_icu_days_unadj, c("peth_binned=High")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with hospital days as an outcome
fit_hosp_days_unadj <- with(imp5, rms::orm(total_number_of_days_in_ho ~ peth_binned))
tbl_hosp_days_unadj <- my_pool_summary(fit_hosp_days_unadj, c("peth_binned=High")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with in-hospital mortality as an outcome
fit_mortality_unadj <- with(imp5, rms::orm(patient_discharged_to_morge ~ peth_binned))
tbl_mortality_unadj <- my_pool_summary(fit_mortality_unadj, c("peth_binned=High")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with discharge as an outcome
fit_discharge_unadj <- with(imp5, rms::orm(patient_discharged_to_home ~ peth_binned))
tbl_discharge_unadj <- my_pool_summary(fit_discharge_unadj, c("peth_binned=High")) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

##### Fit adjusted ORMs to multiply imputed datasets #####

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
tbl_mv_days_adj <- my_pool_summary(fit_mv_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with ICU days as an outcome
fit_icu_days_adj <- with(
  imp5,
  rms::orm(total_no_icu_days ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_icu_days_adj <- my_pool_summary(fit_icu_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with hospital days as an outcome
fit_hosp_days_adj <- with(
  imp5,
  rms::orm(total_number_of_days_in_ho ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_hosp_days_adj <- my_pool_summary(fit_hosp_days_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with in-hospital mortality as an outcome
fit_mortality_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_morge ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_mortality_adj <- my_pool_summary(fit_mortality_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

# perform pooled analysis with discharge as an outcome
fit_discharge_adj <- with(
  imp5,
  rms::orm(patient_discharged_to_home ~ age + sex_male + race_nonwhite + hispanic +
             apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker +
             covid_cat + peth_binned))
tbl_discharge_adj <- my_pool_summary(fit_discharge_adj, eff_names, eff_labels) %>%
  mutate(`odds ratio`=exp(coef), `95% CI lwr`=exp(ci_lwr), `95% CI upr`=exp(ci_upr)) %>%
  mutate(`p-value`=pvalr(p.value)) %>%
  select(effect, `odds ratio`, `95% CI lwr`, `95% CI upr`, `p-value`)

##### Table 3: Effects of alcohol misuse on selected outcomes (adj/unadj) #####

# begin tabulating unadjusted effects of PEth
tbl_mv_days_unadj$Outcome <- "Mechanical ventilation days"
tbl_icu_days_unadj$Outcome <- "Intensive care unit days"
tbl_hosp_days_unadj$Outcome <- "Hospital days"
tbl_mortality_unadj$Outcome <- "Mortality"
tbl_discharge_unadj$Outcome <- "Discharge to home"
tbl_all_unadj <- bind_rows(
  tbl_mv_days_unadj,
  tbl_icu_days_unadj,
  tbl_hosp_days_unadj,
  tbl_mortality_unadj,
  tbl_discharge_unadj) %>%
  select(Outcome, `odds ratio`, contains("95%"), `p-value`)

# begin tabulating adjusted effects of PEth
tbl_mv_days_adj$Outcome <- "Mechanical ventilation days"
tbl_icu_days_adj$Outcome <- "Intensive care unit days"
tbl_hosp_days_adj$Outcome <- "Hospital days"
tbl_mortality_adj$Outcome <- "Mortality"
tbl_discharge_adj$Outcome <- "Discharge to home"
tbl_all_adj <- bind_rows(
  filter(tbl_mv_days_adj, effect=="PEth level: High"),
  filter(tbl_icu_days_adj, effect=="PEth level: High"),
  filter(tbl_hosp_days_adj, effect=="PEth level: High"),
  filter(tbl_mortality_adj, effect=="PEth level: High"),
  filter(tbl_discharge_adj, effect=="PEth level: High")) %>%
  select(Outcome, `odds ratio`, contains("95%"), `p-value`)

# round ORs and CIs to two decimals
# tbl_all_unadj[, 2:4] <- round(tbl_all_unadj[, 2:4], 2)
# tbl_all_adj[, 2:4] <- round(tbl_all_adj[, 2:4], 2)

## EDITED BY RYAN 2023-11-08

##### CC sensitivity analysis, Fit unadjusted ORMs #####

# perform pooled analysis with days on MV as an outcome
my_orm_summ <- function(fit, eff_name) {
  estimates <- fit$coefficients[eff_name]
  vars <- fit$var[eff_name, eff_name]
  n <- fit$stats["Obs"]
  z_stat <- estimates/sqrt(vars)
  p_val <- pnorm(abs(z_stat), lower.tail=F) * 2
  return(data.frame(coef=estimates, std.err=sqrt(vars), 
             z.stat=z_stat, p.value=p_val))
}

fit_mv_days_unadj <- with(df_unimputed, rms::orm(number_mv_days ~ peth_binned))
fit_icu_days_unadj <- with(df_unimputed, rms::orm(total_no_icu_days ~ peth_binned))
fit_hosp_days_unadj <- with(df_unimputed, rms::orm(total_number_of_days_in_ho ~ peth_binned))
fit_mortality_unadj <- with(df_unimputed, rms::orm(patient_discharged_to_morge ~ peth_binned))
fit_discharge_unadj <- with(df_unimputed, rms::orm(patient_discharged_to_home ~ peth_binned))

fit_unadj_cc <- rbind(
  "MV Days" = my_orm_summ(fit_mv_days_unadj, c("peth_binned=High")), 
  "ICU Days" = my_orm_summ(fit_icu_days_unadj, c("peth_binned=High")),
  "Hospital Days" = my_orm_summ(fit_hosp_days_unadj, c("peth_binned=High")),
  "Mortality" = my_orm_summ(fit_mortality_unadj, c("peth_binned=High")),
  "Discharge to home" = my_orm_summ(fit_discharge_unadj, c("peth_binned=High"))
)

fit_mv_days_adj <- with(df_unimputed, rms::orm(number_mv_days ~ age + sex_male + race_nonwhite + hispanic + apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker + covid_cat + peth_binned))
fit_icu_days_adj <- with(df_unimputed, rms::orm(total_no_icu_days ~ age + sex_male + race_nonwhite + hispanic + apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker + covid_cat +peth_binned))
fit_hosp_days_adj <- with(df_unimputed, rms::orm(total_number_of_days_in_ho ~ age + sex_male + race_nonwhite + hispanic + apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker + covid_cat +peth_binned))
fit_mortality_adj <- with(df_unimputed, rms::orm(patient_discharged_to_morge ~ age + sex_male + race_nonwhite + hispanic + apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker + covid_cat +peth_binned))
fit_discharge_adj <- with(df_unimputed, rms::orm(patient_discharged_to_home ~ age + sex_male + race_nonwhite + hispanic + apache_2_admit + diabetes_y_n + cirrhosis_y_n + current_smoker + covid_cat +peth_binned))

fit_adj_cc <- rbind(
  "MV Days" = my_orm_summ(fit_mv_days_adj, c("peth_binned=High")), 
  "ICU Days" = my_orm_summ(fit_icu_days_adj, c("peth_binned=High")),
  "Hospital Days" = my_orm_summ(fit_hosp_days_adj, c("peth_binned=High")),
  "Mortality" = my_orm_summ(fit_mortality_adj, c("peth_binned=High")),
  "Discharge to home" = my_orm_summ(fit_discharge_adj, c("peth_binned=High"))
)
