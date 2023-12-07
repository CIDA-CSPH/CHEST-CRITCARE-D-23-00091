################################################### -
## Title: Prepare data for analysis
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## PI: Ellen Burnham
## Date Created: 2023-10-19
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(janitor, warn.conflicts=FALSE)
library(readxl)
library(CIDAtools)

# read raw dataset containing demographics, clinical vars, and AUD vars
df_raw <- read.csv(ProjectLocation("DataRaw/data-2022-07-29.csv")) %>%
  clean_names() %>%
  filter(!grepl("\\*", record_id)) %>% ### removes patients who did not consent
  filter(record_id!="") ###removes blank rows at end of file

# read dataset of corrected race/ethnicity
df_revisions <- read_xlsx(ProjectLocation("DataRaw/RF_missing race_090222.xlsx")) %>%
  clean_names() %>%
  rename(revised_hispanic=hispanic)

# read dataset of additional PEth measures
df_peth_addl <- read_xlsx(
  ProjectLocation("DataRaw/Day 1 RF RBC's that were previously missing.xlsx"),
  range="A2:B17") %>%
  rename(revised_peth=`Amount *(ng/mL)`) %>%
  mutate(record_id=gsub(" Day 1", "", PID))

# read dataset of delirium variables
df_delirium <- read.csv(ProjectLocation("DataRaw/Delirium_RF_retro_100422.csv")) %>%
  clean_names()

# NOTE: subjects with delirium data not represented in main dataset
ids_extra <- df_delirium$record_id[!(df_delirium$record_id %in% df_raw$record_id)]
warning(paste0(length(ids_extra), " subjects with delirium are not in main file..."))

##### Merge data sources #####

# merge new race and ethnicity labels, using newer version if available
df_merged <- df_raw %>%
  left_join(select(df_revisions, record_id, revised_hispanic, revised_race_cat)) %>%
  mutate(
    hispanic=coalesce(revised_hispanic, hispanic),
    race=coalesce(revised_race_cat, race))

# merge new PEth measures, using newer version if available
df_merged <- df_merged %>%
  left_join(select(df_peth_addl, record_id, revised_peth)) %>%
  mutate(peth=coalesce(peth, revised_peth))

# merge delirium data
df_merged <- df_merged %>%
  left_join(select(df_delirium, record_id, contains("_delirium"))) %>%
  select(-one_date_that_delirium_was)

##### Derive additional variables #####

# classify pre- / peri- pandemic by year of enrollment
df_merged <- df_merged %>%
  mutate(
    pandemic_phase=ifelse(year_enrolled <= 2019, "Pre-pandemic", "Peri-pandemic"),
    prepandemic=ifelse(year_enrolled <= 2019, 1, 0),
    peripandemic=ifelse(year_enrolled > 2019, 1, 0)) %>%
  mutate(pandemic_phase=factor(pandemic_phase, c("Pre-pandemic", "Peri-pandemic"))) %>%
  mutate(covid_cat=factor(case_when(
    covid=="Y" ~ "Positive",
    covid=="N" & peripandemic ~ "Negative",
    covid=="N" & prepandemic ~ "Pre-pandemic"),
    c("Negative", "Positive", "Pre-pandemic")))

# categorize demographics as in previous analysis
df_merged <- df_merged %>%
  mutate(sex_male=ifelse(gender==1, TRUE, FALSE)) %>%
  mutate(hispanic=case_when(
    hispanic==1 ~ TRUE,
    hispanic==0 ~ FALSE)) %>%
  mutate(race_nonwhite=case_when(
    race==5 ~ FALSE,
    race %in% c(1, 2, 3, 4, 6) ~ TRUE)) %>%
  mutate(race=case_when(
    race==1 ~ "American Indian",
    race==2 ~ "Asian",
    race==3 ~ "African American",
    race==4 ~ "Pacific Islander",
    race==5 ~ "White",
    race==6 ~ "More than one race",
    TRUE ~ "Unknown / Refused to answer")) %>%
  mutate(diabetes_y_n=case_when(
    diabetes_y_n=="Checked" ~ TRUE,
    diabetes_y_n=="Unchecked" ~ FALSE)) %>%
  mutate(cirrhosis_y_n=case_when(
    cirrhosis_y_n=="Checked" ~ TRUE,
    cirrhosis_y_n=="Unchecked" ~ FALSE)) %>%
  mutate(current_smoker=case_when(
    current_smoker=="Yes" ~ TRUE,
    current_smoker=="No" ~ FALSE)) %>%
  mutate(likely_current_drinking=case_when(
    likely_current_drinking==1 ~ TRUE,
    likely_current_drinking==0 ~ FALSE)) %>%
  mutate(any_aud_all_data=case_when(
    any_aud_all_data==0 ~ FALSE,
    any_aud_all_data==1 ~ TRUE,
    any_aud_all_data==2 ~ NA))

# categorize smoking history
df_merged <- df_merged %>%
  mutate(smoking_history=case_when(
    current_smoker ~ "Current smoker",
    !current_smoker ~ "Never or former smoker"))

# bin PEth measures according to AFSHAR ET AL. (2022)
df_merged <- df_merged %>%
  mutate(peth_binned = factor(case_when(
    peth < 25 ~ "Low",
    peth >= 25 ~ "High"),
    c("Low", "High")))

# bin PEth measures according to ULWELLING ET AL. (2018)
df_merged <- df_merged %>%
  mutate(peth_binned_ulwelling=factor(case_when(
    peth < 20 ~ "Light/none",
    peth < 200 ~ "Significant",
    peth >= 200 ~ "Heavy"),
    c("Light/none", "Significant", "Heavy")))

# compute log-transformed PEth
df_merged$peth_logged <- log(df_merged$peth + 10)

# derive AUD based on RN-validated AUDIT-C
df_merged <- df_merged %>%
  mutate(any_aud_audit_c=case_when(
    audit_c_nurse_validated >= 4 & sex_male ~ 1,
    audit_c_nurse_validated >= 3 & !sex_male ~ 1,
    audit_c_nurse_validated < 5 ~ 0))

# bin delirium categories
df_merged <- df_merged %>%
  mutate(delirium_grouped = factor(case_when(
    definitive_delirium==1 | probable_delirium==1 |
      possible_delirium==1 ~ "Possible/Probable/Definitive",
    no_delirium==1 ~ "None"))) 

# impute one missing covid case in the prepandemic phase
df_merged <- df_merged %>%
  mutate(covid=ifelse(covid=="", "N", covid))

##### Write results of data preparation to output file #####

# define output file name
out_name <- ProjectLocation("DataProcessed/merged_analysis_dataset.csv")

# write to output file
write.csv(df_merged, out_name, row.names=F)
message("File written to: ", out_name)
