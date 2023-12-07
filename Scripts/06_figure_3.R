################################################### -
## Title: Analysis of PEth and AUDIT-C (for CHEST revision)
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## PI: Ellen Burnham
## Date Created: 2023-11-01
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(mice, exclude="filter")
library(miceadds)
library(ggplot2)
theme_set(ggpubr::theme_classic2(base_size=15))
library(CIDAtools)
library(here)

# load multiply-imputed data
imp5 <- readRDS(ProjectLocation("DataProcessed/multiply_imputed_dataset.rds"))

# isolate original (unimputed) data
df_unimputed <- complete(imp5, 0)

# plot Log-Scale PEth versus AUDIT-C (with pandemic interaction)
fig_peth_audit_strat <- df_unimputed %>%
  filter(!is.na(audit_c_nurse_validated) & !is.na(peth_logged)) %>%
  ggplot(aes(audit_c_nurse_validated, peth_logged, col=pandemic_phase)) +
  geom_smooth(method="lm", formula="y ~ x", se=FALSE) +
  geom_jitter(width=0.2, height=0.2, alpha=0.5) +
  labs(x="AUDIT-C (RN validated)", y="Log-scale PEth", col="Phase",
       title="Log-scale PEth vs. AUDIT-C")

# complete case analysis
cor.test(
  df_unimputed$peth_logged[df_unimputed$pandemic_phase=="Pre-pandemic"],
  df_unimputed$audit_c_nurse_validated[df_unimputed$pandemic_phase=="Pre-pandemic"],
  use="complete.obs")
cor.test(
  df_unimputed$peth_logged[df_unimputed$pandemic_phase=="Peri-pandemic"],
  df_unimputed$audit_c_nurse_validated[df_unimputed$pandemic_phase=="Peri-pandemic"],
  use="complete.obs")

# compute correlation between PEth and AUDIT-C
imp5_prepandemic <- subset_datlist(imp5, imp5[[1]]$pandemic_phase=="Pre-pandemic")
cor_prepandemic <- micombine.cor(mi.res=imp5_prepandemic, variables=c(15, 30))
r_prepandemic <- round(cor_prepandemic$r[1], 2)
p_prepandemic <- format.pval(cor_prepandemic$p[1], digits=2, eps=0.01)

imp5_peripandemic <- subset_datlist(imp5, imp5[[1]]$pandemic_phase=="Peri-pandemic")
cor_peripandemic <- micombine.cor(mi.res=imp5_peripandemic, variables=c(15, 30))
r_peripandemic <- round(cor_peripandemic$r[1], 2)
p_peripandemic <- format.pval(cor_peripandemic$p[1], digits=2, eps=0.01)

# annotate figure with correlations
fig_peth_audit_strat <- fig_peth_audit_strat +
  annotate("text", x=8, y=3, col="#F8766D", label=paste0(
    "Pearson Rho = ", r_prepandemic, ", p ", p_prepandemic)) +
  annotate("text", x=4, y=9, col="#00BFC4", label=paste0(
    "Pearson Rho = ", r_peripandemic, ", p ", p_peripandemic))

ggsave(plot=fig_peth_audit_strat, filename=here("Dissemination/Figure3.png"), width=8, height=6)