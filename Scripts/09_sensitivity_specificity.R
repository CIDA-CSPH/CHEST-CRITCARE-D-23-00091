################################################### -
## Title: Analysis of PEth vs. AUDIT-C Sensitivity / Specificity
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## PI: Ellen Burnham
## Date Created: 2023-11-16
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(mice, exclude="filter")
library(CIDAtools)

# load multiply-imputed data
imp5 <- readRDS(ProjectLocation("DataProcessed/multiply_imputed_dataset.rds"))

# isolate original (unimputed) data
df_unimputed <- complete(imp5, 0)

##### Function to calculate operating characteristics #####

my.performance <- function(df) {
  mat_conf <- table(df$any_aud_audit_c, df$peth_binned)
  a <- mat_conf[2, 2]
  b <- mat_conf[2, 1]
  c <- mat_conf[1, 2]
  d <- mat_conf[1, 1]
  sens <- a / (a + c)
  spec <- d / (b + d)
  prev <- (a + c)/(a + b + c + d)
  ppv <- (sens * prev)/((sens * prev) + ((1 - spec) * (1 - prev)))
  npv <- (spec * (1 - prev))/(((1- sens) * prev) + ((spec) * (1 - prev)))
  list(sensitivity=sens, specificity=spec, PPV=ppv, NPV=npv)
}

##### Function to pool estimates across multiply imputed datasets #####

# use ML theory and Rubin's rules to estimate CI for sensitivity/specificty/ppv/npv
my.pooler <- function(df, res, col=1){
  n <- nrow(df)
  pooled <- pool.scalar(
    Q=res[, col],
    U=(res[, col] * (1 - res[, col])) / n,
    n=n, k=1)
  stat <- pooled$qbar / sqrt(pooled$t)
  ci_95 <- pooled$qbar + sqrt(pooled$t) * c(-1.96, 1.96)
  list(mean=pooled$qbar, lwr=ci_95[1], upr=ci_95[2])
}

##### Calculate operating characteristics across pandemic eras #####

# start with point esitimates (both pre-pandemic and post-pandemic)
df_prepandemic <- filter(df_unimputed, pandemic_phase=="Pre-pandemic")
df_peripandemic <- filter(df_unimputed, pandemic_phase=="Peri-pandemic")

# store results
res_prepandemic <- data.frame(my.performance(df_prepandemic))
res_peripandemic <- data.frame(my.performance(df_peripandemic))

# iterate over multiply imputed datasets
for (i in 1:5) {
  df_prepandemic <- filter(complete(imp5, i), pandemic_phase=="Pre-pandemic")
  df_peripandemic <- filter(complete(imp5, i), pandemic_phase=="Peri-pandemic")
  
  tmp_prepandemic <- data.frame(my.performance(df_prepandemic))
  tmp_peripandemic <- data.frame(my.performance(df_peripandemic))
  
  res_prepandemic <- bind_rows(res_prepandemic, tmp_prepandemic)
  res_peripandemic <- bind_rows(res_peripandemic, tmp_peripandemic)
}

# remove results from unimputed datasets (only used for sanity check)
res_prepandemic <- res_prepandemic[-1, ]
res_peripandemic <- res_peripandemic[-1, ]

# present in tabular format
means_pre <- rep(NA, 4)
cis_pre <- rep(NA, 4)
means_peri <- rep(NA, 4)
cis_peri <- rep(NA, 4)
for (i in 1:4) {
  pool_pre <- my.pooler(df_prepandemic, res_prepandemic, i)
  pool_peri <- my.pooler(df_peripandemic, res_peripandemic, i)
  
  means_pre[i] <- round(pool_pre$mean, 3)
  means_peri[i] <- round(pool_peri$mean, 3)
  cis_pre[i] <- paste0("(", round(pool_pre$lwr, 2), ", ", round(pool_pre$upr, 2), ")")
  cis_peri[i] <- paste0("(", round(pool_peri$lwr, 2), ", ", round(pool_peri$upr, 2), ")")
}

tab_sens_spec <- data.frame(
  Metric=c("Sensitivity", "Specificity", "PPV", "NPV"),
  Prepandemic=means_pre,
  `95% Conf. Int`=cis_pre,
  Peripandemic=means_peri,
  `95% Conf. Int`=cis_peri,
  check.names=FALSE)


