################################################### -
## Title: Create figures for supplementary materials
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## PI: Ellen Burnham
## Date Created: 2023-11-01
################################################### -


options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw(base_size=10))
library(CIDAtools)
library(here)

# read analysis dataset
df_analysis <- read.csv(ProjectLocation("DataProcessed/merged_analysis_dataset.csv"))

##### Supplemental Figure 1: AUDIT-C scores among critically ill patients over time #####

# Log-scale PEth measures over time
fig_s1 <- df_analysis %>%
  mutate(sex=ifelse(gender==1, "Male", "Female")) %>%
  ggplot(aes(year_enrolled, audit_c_nurse_validated, color=sex)) +
  geom_jitter(width=0.1, height=0, alpha=0.45, na.rm=T) +
  geom_smooth(se=F, na.rm=T, method="loess", formula="y ~ x") +
  geom_vline(xintercept=2020, linetype="dashed") +
  labs(x="Year enrolled", y="AUDIT-C (RN Validated)", color="Sex",
       title="AUDIT scores over time (by Sex)")

ggsave(plot=fig_s1, filename=here("Dissemination/SupplementalFigure1.png"), width=8, height=6)

