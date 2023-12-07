################################################### -
## Title: Create manuscript Figure 2
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

# factor trichotomous PEth levels
df_analysis <- df_analysis %>%
  mutate(peth_binned_ulwelling=factor(
    peth_binned_ulwelling,
    c("Light/none", "Significant", "Heavy"))) %>%
  mutate(`PEth level`=factor(
    peth_binned_ulwelling,
    c("Light/none", "Significant", "Heavy"),
    c("<20 ng/mL", "20-200 ng/mL", ">=200 ng/mL")))

##### Figure 2. Relationship between patient PEth values and demographics #####

fig_age_peth <- df_analysis %>%
  filter(!is.na(peth_binned_ulwelling) & !is.na(age)) %>%
  ggboxplot(
    x="peth_binned_ulwelling",
    y="age",
    fill="PEth level",
    xlab="",
    ylab="Age",
    add=c("dotplot")) +
  stat_compare_means(label.y=95)

fig_bmi_peth <- df_analysis %>%
  filter(!is.na(peth_binned_ulwelling) & !is.na(body_mass_index)) %>%
  ggboxplot(
    x="peth_binned_ulwelling",
    y="body_mass_index",
    fill="PEth level",
    xlab="",
    ylab="BMI",
    add=c("dotplot")) +
  stat_compare_means(label.y=85)

pval_sex_peth <- format.pval(
  fisher.test(
    table(df_analysis$sex_male, df_analysis$peth_binned_ulwelling))$p.value,
  digits=3)
fig_sex_peth <- df_analysis %>%
  filter(!is.na(peth_binned_ulwelling)) %>%
  mutate(sex=ifelse(sex_male==1, "Male", "Female")) %>%
  group_by(peth_binned_ulwelling, sex) %>%
  summarise(count=n()) %>%
  mutate(pct=prop.table(count)*100) %>%
  ggplot(aes(fill=sex, y=count, x=peth_binned_ulwelling)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5), angle=45, size=3) +
  labs(x="PEth level", y="Number of patients", fill="Sex") +
  theme_classic2() +
  annotate("text", x=2.5, y=75, label=paste0("Fisher's Exact, p = ", pval_sex_peth))

pval_hispanic_peth <- format.pval(
  fisher.test(
    table(df_analysis$hispanic, df_analysis$peth_binned_ulwelling))$p.value,
  digits=3)
fig_hispanic_peth <- df_analysis %>%
  filter(!is.na(peth_binned_ulwelling) & !is.na(hispanic)) %>%
  mutate(ethnicity=ifelse(hispanic==1, "Hispanic", "Non-Hispanic")) %>%
  group_by(peth_binned_ulwelling, ethnicity) %>%
  summarise(count=n()) %>%
  mutate(pct=prop.table(count)*100) %>%
  ggplot(aes(fill=ethnicity, y=count, x=peth_binned_ulwelling)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5), angle=45, size=3) +
  labs(x="PEth level", y="Number of patients", fill="Ethnicity") +
  theme_classic2() +
  annotate("text", x=2.5, y=75, label=paste0("Fisher's Exact, p = ", pval_hispanic_peth))

pval_race_peth <- format.pval(
  fisher.test(
    table(df_analysis$race_nonwhite, df_analysis$peth_binned_ulwelling))$p.value,
  digits=3)
fig_race_peth <- df_analysis %>%
  filter(!is.na(peth_binned_ulwelling) & !is.na(race_nonwhite)) %>%
  mutate(race=ifelse(race_nonwhite==1, "Non-White", "White")) %>%
  group_by(peth_binned_ulwelling, race) %>%
  summarise(count=n()) %>%
  mutate(pct=prop.table(count)*100) %>%
  ggplot(aes(fill=race, y=count, x=peth_binned_ulwelling)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5), angle=45, size=3) +
  labs(x="PEth level", y="Number of patients", fill="Race") +
  theme_classic2() +
  annotate("text", x=2.5, y=75, label=paste0("Fisher's Exact, p = ", pval_race_peth))

fig_2 <- ggarrange(
  fig_age_peth, fig_bmi_peth, "",
  fig_sex_peth, fig_hispanic_peth, fig_race_peth,
  labels=c("A", "B", "", "C", "D", "E"),
  ncol=3, nrow=2)

ggsave(plot=fig_2, filename=here("Dissemination/Figure2.png"), width=14, height=10)
