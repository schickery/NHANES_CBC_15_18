install.packages("dplyr")
install.packages("tidyr")
install.packages("survey")
install.packages("ggplot2")
install.packages("gridExtra")
library(dplyr)
library(tidyr)
library(survey)
library(ggplot2)
library(gridExtra)

#Install Files from NHANES Database
#Demographics files
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/DEMO_H.XPT", tf <- tempfile(), mode="wb")
DEMO_H <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDEXPRG", "RIDRETH3", "SDMVSTRA","SDMVPSU","WTMEC2YR")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT", tf <- tempfile(), mode="wb")
DEMO_I <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDEXPRG", "RIDRETH3", "SDMVSTRA","SDMVPSU","WTMEC2YR")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_DEMO.XPT", tf <- tempfile(), mode="wb")
DEMO_P <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDEXPRG", "RIDRETH3", "SDMVSTRA","SDMVPSU","WTMECPRP")]


#CBC Data (CBC)
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/CBC_H.XPT", tf <- tempfile(), mode="wb")
CBC_H <- foreign::read.xport(tf)[,c("SEQN", "LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/CBC_I.XPT", tf <- tempfile(), mode="wb")
CBC_I <- foreign::read.xport(tf)[,c("SEQN", "LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_CBC.XPT", tf <- tempfile(), mode="wb")
CBC_P <- foreign::read.xport(tf)[,c("SEQN", "LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")]

#Biochem for Creatinine and Glucose
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BIOPRO_H.XPT", tf <- tempfile(), mode="wb")
BIOPRO_H <- foreign::read.xport(tf)[,c("SEQN", "LBXSCR", "LBXSGL")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BIOPRO_I.XPT", tf <- tempfile(), mode="wb")
BIOPRO_I <- foreign::read.xport(tf)[,c("SEQN", "LBXSCR", "LBXSGL")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_BIOPRO.XPT", tf <- tempfile(), mode="wb")
BIOPRO_P <- foreign::read.xport(tf)[,c("SEQN", "LBXSCR", "LBXSGL")]

#Hepatitis B Surface Ag and Core Ab
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/HEPBD_H.XPT", tf <- tempfile(), mode="wb")
HEPBD_H <- foreign::read.xport(tf)[,c("SEQN", "LBXHBC", "LBDHBG")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HEPBD_I.XPT", tf <- tempfile(), mode="wb")
HEPBD_I <- foreign::read.xport(tf)[,c("SEQN", "LBXHBC", "LBDHBG")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_HEPBD.XPT", tf <- tempfile(), mode="wb")
HEPBD_P <- foreign::read.xport(tf)[,c("SEQN", "LBXHBC", "LBDHBG")]

#HEP C Confirmed Data
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/HEPC_H.XPT", tf <- tempfile(), mode="wb")
HEPC_H <- foreign::read.xport(tf)[,c("SEQN", "LBXHCR")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_HEPC.XPT", tf <- tempfile(), mode="wb")
HEPC_P <- foreign::read.xport(tf)[,c("SEQN", "LBXHCR")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HEPC_I.XPT", tf <- tempfile(), mode="wb")
HEPC_I <- foreign::read.xport(tf)[,c("SEQN", "LBXHCR")]

#Blood Donation & Illness (Exclude Past 12 weeks donators)
#download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/HSQ_H.XPT", tf <- tempfile(), mode="wb")
#HSQ_H <- foreign::read.xport(tf)[,c("SEQN", "HSQ500", "HSQ510", "HSQ520","HSQ571", "HSQ580")]
#download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSQ_I.XPT", tf <- tempfile(), mode="wb")
#HSQ_I <- foreign::read.xport(tf)[,c("SEQN", "HSQ500", "HSQ510", "HSQ520","HSQ571", "HSQ580")]


#Alcohol Use
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/ALQ_H.XPT", tf <- tempfile(), mode="wb")
ALQ_H <- foreign::read.xport(tf)[,c("SEQN", "ALQ141Q", "ALQ141U")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/ALQ_I.XPT", tf <- tempfile(), mode="wb")
ALQ_I <- foreign::read.xport(tf)[,c("SEQN", "ALQ141Q", "ALQ141U")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_ALQ.XPT", tf <- tempfile(), mode="wb")
ALQ_P <- foreign::read.xport(tf)[,c("SEQN", "ALQ142")]

#Smoking
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/SMQ_H.XPT", tf <- tempfile(), mode="wb")
SMQ_H <- foreign::read.xport(tf)[,c("SEQN", "SMQ040")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SMQ_I.XPT", tf <- tempfile(), mode="wb")
SMQ_I <- foreign::read.xport(tf)[,c("SEQN", "SMQ040")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_SMQ.XPT", tf <- tempfile(), mode="wb")
SMQ_P <- foreign::read.xport(tf)[,c("SEQN", "SMQ040")]

#A1C
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/GHB_H.XPT", tf <- tempfile(), mode="wb")
GHB_H <- foreign::read.xport(tf)[,c("SEQN", "LBXGH")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/GHB_I.XPT", tf <- tempfile(), mode="wb")
GHB_I <- foreign::read.xport(tf)[,c("SEQN", "LBXGH")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_GHB.XPT", tf <- tempfile(), mode="wb")
GHB_P <- foreign::read.xport(tf)[,c("SEQN", "LBXGH")]

#hsCRP >7.5 MG/L
#download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSCRP_I.XPT", tf <- tempfile(), mode="wb")
#HSCRP_I <- foreign::read.xport(tf)[,c("SEQN", "LBXHSCRP")]
#download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_HSCRP.XPT", tf <- tempfile(), mode="wb")
#HSCRP_P <- foreign::read.xport(tf)[,c("SEQN", "LBXHSCRP")]

#Pregnancy and Breast Feeding
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/RHQ_H.XPT", tf <- tempfile(), mode="wb")
RHQ_H <- foreign::read.xport(tf)[,c("SEQN", "RHD143", "RHQ200")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/RHQ_I.XPT", tf <- tempfile(), mode="wb")
RHQ_I <- foreign::read.xport(tf)[,c("SEQN", "RHD143", "RHQ200")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_RHQ.XPT", tf <- tempfile(), mode="wb")
RHQ_P <- foreign::read.xport(tf)[,c("SEQN", "RHD143", "RHQ200")]

#BMI greater than 30 BMXBMI
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BMX_H.XPT", tf <- tempfile(), mode="wb")
BMX_H <- foreign::read.xport(tf)[,c("SEQN", "BMXBMI")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BMX_I.XPT", tf <- tempfile(), mode="wb")
BMX_I <- foreign::read.xport(tf)[,c("SEQN", "BMXBMI")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_BMX.XPT", tf <- tempfile(), mode="wb")
BMX_P <- foreign::read.xport(tf)[,c("SEQN", "BMXBMI")]

#Append Files
DEMO <- bind_rows(DEMO_H, DEMO_I, DEMO_P)
CBC <- bind_rows(CBC_H, CBC_I, CBC_P)
BIOPRO <- bind_rows(BIOPRO_H, BIOPRO_I, BIOPRO_P)
HEPBD <- bind_rows(HEPBD_H, HEPBD_I, HEPBD_P)
HEPC <- bind_rows(HEPC_H, HEPC_I, HEPC_P)
#HSQ <- bind_rows(HSQ_H, HSQ_I)
ALQ <- bind_rows(ALQ_H, ALQ_I, ALQ_P)
SMQ <- bind_rows(SMQ_H, SMQ_I, SMQ_P)
GHB <- bind_rows(GHB_H, GHB_I, GHB_P)
RHQ <- bind_rows(RHQ_H, RHQ_I, RHQ_P)
BMX <- bind_rows(BMX_H, BMX_I, BMX_P)

#------------------------------------------
#Determine number of SEQNs before exclusions

install.packages("purrr")
library(purrr)

# Combine all datasets into one dataframe
combined_data <- reduce(list(DEMO, CBC, BIOPRO, HEPBD, HEPC, ALQ, SMQ, GHB, RHQ, BMX), full_join, by = "SEQN")

# Check the number of rows in the combined dataset
number_of_rows <- nrow(combined_data)
print(number_of_rows)

# Adjust weights for combined periods
combined_data <- combined_data %>%
  mutate(WTMEC4YR = case_when(
    between(SEQN, min(DEMO_H$SEQN), max(DEMO_H$SEQN)) ~ WTMEC2YR * (2 / 7.2),
    between(SEQN, min(DEMO_I$SEQN), max(DEMO_I$SEQN)) ~ WTMEC2YR * (2 / 7.2),
    between(SEQN, min(DEMO_P$SEQN), max(DEMO_P$SEQN)) ~ WTMECPRP * (3.2 / 7.2),
    TRUE ~ NA_real_
  ))

# Check the first few rows to ensure the weights are adjusted
head(combined_data)

# Create a survey design object using the combined data with adjusted weights
nhanes_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC4YR,
  data = combined_data,
  nest = TRUE
)

# Check the summary of the survey design object
summary(nhanes_design)

# Count the number of unique SEQN values
unique_seqn_count <- n_distinct(combined_data$SEQN)
print(paste("Number of unique SEQN in the combined dataset:", unique_seqn_count))

summarize_counts <- function(data, design, file_name) {
  # Summarizing counts by gender and race using the survey design
  gender_race_weighted_counts <- svytable(~RIAGENDR + RIDRETH3, design)
  
  # Convert to a data frame
  gender_race_weighted_df <- as.data.frame(gender_race_weighted_counts)
  colnames(gender_race_weighted_df) <- c("gender", "race", "count")
  
  # Decode gender and race
  gender_race_weighted_df <- gender_race_weighted_df %>%
    mutate(
      gender = recode(gender, `1` = "Male", `2` = "Female"),
      race = recode(race,
                    `1` = "Mexican American",
                    `2` = "Other Hispanic",
                    `3` = "Non-Hispanic White",
                    `4` = "Non-Hispanic Black",
                    `6` = "Non-Hispanic Asian",
                    `7` = "Other Race - Including Multi-Racial")
    )
  
  # Pivot the data to wide format
  gender_race_table <- gender_race_weighted_df %>%
    pivot_wider(names_from = race, values_from = count, values_fill = list(count = 0))
  
  # Print the table
  print(gender_race_table)
  
  # Create a table grob
  table_grob <- tableGrob(gender_race_table)
  
  # Adjust column widths
  col_widths <- unit(rep(1, ncol(table_grob)), "npc") / ncol(table_grob)
  table_grob$widths <- col_widths
  
  # Save the table grob as a PDF with adjusted page size
  pdf(file_name, width = 12, height = 6)  # Adjust width and height as needed
  grid.draw(table_grob)
  dev.off()
}

# Initial summary without exclusions
summarize_counts(combined_data, nhanes_design, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/NHANES_CBC_15_18/initial_summary.pdf")

# Apply exclusion criteria: A1C >= 6.5%--------------------------------------------------------------------------------------
combined_data_excluded_a1c <- combined_data %>%
  filter(LBXGH < 6.5 | is.na(LBXGH))

# Create a new survey design object using the filtered data
nhanes_design_excluded_a1c <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_a1c, nhanes_design_excluded_a1c, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion.pdf")

# Apply exclusion criteria: Serum creatinine > 2.5 mg/dL--------------------------------------------------
combined_data_excluded_creatinine <- combined_data_excluded_a1c %>%
  filter(LBXSCR < 2.5 | is.na(LBXSCR))

# Create a new survey design object using the filtered data
nhanes_design_excluded_creatinine <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_creatinine, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_creatinine, nhanes_design_excluded_creatinine, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_creatinine.pdf")

# Apply exclusion criteria: Hepatitis C positivity (LBXHCR == 1)--------------------------------------------------------
combined_data_excluded_hepc <- combined_data_excluded_creatinine %>%
  filter(LBXHCR != 1 | is.na(LBXHCR))

# Create a new survey design object using the filtered data
nhanes_design_excluded_hepc <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_hepc, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_hepc, nhanes_design_excluded_hepc, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_hepc.pdf")

# Apply exclusion criteria: Hepatitis B surface antigen (LBDHBG == 1) or core antibody (LBXHBC == 1)---------------------------------
combined_data_excluded_hepb <- combined_data_excluded_hepc %>%
  filter((LBDHBG != 1 | is.na(LBDHBG)) & (LBXHBC != 1 | is.na(LBXHBC)))

# Create a new survey design object using the filtered data
nhanes_design_excluded_hepb <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_hepb, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_hepb, nhanes_design_excluded_hepb, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_hepb.pdf")

# Apply exclusion criteria: Blood donation in the past 12 weeks (HSQ580 == 1, 2, or 3)
#combined_data_excluded_blood_donation <- combined_data_excluded_hepb %>%
#  filter(!(HSQ580 %in% c(1, 2, 3)))

# Create a new survey design object using the filtered data
#nhanes_design_excluded_blood_donation <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_blood_donation, nest = TRUE)

# Summarize counts after exclusion and save as PDF
#summarize_counts(combined_data_excluded_blood_donation, nhanes_design_excluded_blood_donation, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_blood_donation.pdf")

# Apply exclusion criteria for 2015-2016 data: ALQ141Q >= 54----------------------------------------------
combined_data_excluded_alcohol_2015_2016 <- combined_data_excluded_hepb %>%
  filter((is.na(ALQ141Q) | ALQ141Q < 104))

# Apply exclusion criteria for 2017-2018 data: ALQ142 in 1, 2, 3, or 4
combined_data_excluded_alcohol <- combined_data_excluded_alcohol_2015_2016 %>%
  filter(!(ALQ142 %in% c(1, 2, 3, 4)))

# Create a new survey design object using the filtered data
nhanes_design_excluded_alcohol <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_alcohol, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_alcohol, nhanes_design_excluded_alcohol, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_alcohol.pdf")

# Apply exclusion criteria: hsCRP LBXHSCRP >7.5 mg/L
#combined_data_excluded_hscrp <- combined_data_excluded_alcohol %>%
#  filter(LBXHSCRP > 7.5 | is.na(LBXHSCRP))

# Create a new survey design object using the filtered data
#nhanes_design_excluded_hscrp <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_pregnancy_breastfeeding, nest = TRUE)

# Summarize counts after exclusion and save as PDF
#summarize_counts(combined_data_excluded_hscrp, nhanes_design_excluded_hscrp, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_hscrp.pdf")

# Apply exclusion criteria: Pregnancy (RHD143 == 1) and Breastfeeding (RHQ200 == 1) for females
combined_data_excluded_pregnancy_breastfeeding <- combined_data_excluded_alcohol %>%
  filter(!(RIAGENDR == 2 & (RHD143 == 1 | RHQ200 == 1)))

# Create a new survey design object using the filtered data
nhanes_design_excluded_pregnancy_breastfeeding <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_pregnancy_breastfeeding, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_pregnancy_breastfeeding, nhanes_design_excluded_pregnancy_breastfeeding, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_pregnancy_breastfeeding.pdf")

# Apply exclusion criteria: BMI (BMXBMI > 30)
combined_data_excluded_bmi <- combined_data_excluded_pregnancy_breastfeeding %>%
filter(BMXBMI > 30 | is.na(BMXBMI))

# Create a new survey design object using filtered data
nhanes_design_excluded_bmi <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_bmi, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_bmi, nhanes_design_excluded_bmi, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_bmi.pdf")

# Apply exclusion criteria: Age < 20
combined_data_excluded_age <- combined_data_excluded_bmi %>%
  filter(RIDAGEYR > 20)

# Create a new survey design object using the filtered data
nhanes_design_excluded_age <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_age, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_age, nhanes_design_excluded_age, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_age.pdf")

# Remove Missing CBCs
# Define the list of CBC variables
cbc_vars <- c("LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", 
              "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", 
              "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")

# Remove cases with any missing CBC variables
combined_data_complete_cases <- combined_data_excluded_age %>%
  filter(!is.na(LBXWBCSI) & !is.na(LBXLYPCT) & !is.na(LBXMOPCT) & !is.na(LBXNEPCT) & 
           !is.na(LBXEOPCT) & !is.na(LBXBAPCT) & !is.na(LBXRBCSI) & !is.na(LBXHGB) & 
           !is.na(LBXHCT) & !is.na(LBXMCVSI) & !is.na(LBXMCHSI) & !is.na(LBXMC) & 
           !is.na(LBXRDW) & !is.na(LBXPLTSI) & !is.na(LBXMPSI))

# Create a new survey design object using the filtered data
nhanes_design_complete_cases <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_complete_cases, nest = TRUE)

# Check the dimensions and summary of the filtered data
dim(combined_data_complete_cases)
summary(combined_data_complete_cases)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_complete_cases, nhanes_design_complete_cases, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_complete_cases.pdf")

# Define age groups
combined_data_complete_cases <- combined_data_complete_cases %>%
  mutate(age_group = case_when(
    RIDAGEYR >= 20 & RIDAGEYR <= 39 ~ "20-39",
    RIDAGEYR >= 40 & RIDAGEYR <= 59 ~ "40-59",
    RIDAGEYR >= 65 ~ "65+",
    TRUE ~ NA_character_ 
  ))

# Function to filter data by race, gender, and age group
filter_data <- function(data, race, gender, age_group) {
  data %>%
    filter(RIDRETH3 == race & RIAGENDR == gender & age_group == age_group)
}

races <- c(1, 2, 3, 4, 6, 7)  # Replace with actual race codes
genders <- c(1, 2)  # 1 for Male, 2 for Female
age_groups <- c("20-39", "40-59", "65+")

# List to store results
results <- list()

# Loop through each combination
for (race in races) {
  for (gender in genders) {
    for (age_group in age_groups) {
      key <- paste0("race_", race, "_gender_", gender, "_age_", age_group)
      results[[key]] <- filter_data(combined_data_complete_cases, race, gender, age_group)
    }
  }
}

# Count the number of cases for each combination of gender, race, and age group
count_cases <- combined_data_complete_cases %>%
  group_by(RIDRETH3, RIAGENDR, age_group) %>%
  summarise(count = n()) %>%
  ungroup()

# Print the resulting table
print(count_cases)

# Save the summary counts as a CSV file
write.csv(count_cases, file = "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/NHANES_filtered_data_summary.csv", row.names = FALSE)

# Define a function to calculate descriptive statistics
calculate_descriptive_stats <- function(data, variable) {
  data %>%
    summarise(
      N = n(),
      Mean = mean(get(variable), na.rm = TRUE),
      SD = sd(get(variable), na.rm = TRUE),
      Median = median(get(variable), na.rm = TRUE),
      IQR = IQR(get(variable), na.rm = TRUE),
      Min = min(get(variable), na.rm = TRUE),
      Max = max(get(variable), na.rm = TRUE)
    )
}

# Apply this function to the CBC variables grouped by race, gender, and age group
descriptive_stats <- combined_data_complete_cases %>%
  group_by(RIDRETH3, RIAGENDR, age_group) %>%
  summarise(
    across(c("LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", 
             "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", 
             "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI"), 
           list(
             N = ~n(),
             Mean = ~mean(.x, na.rm = TRUE),
             SD = ~sd(.x, na.rm = TRUE),
             Median = ~median(.x, na.rm = TRUE),
             IQR = ~IQR(.x, na.rm = TRUE),
             Min = ~min(.x, na.rm = TRUE),
             Max = ~max(.x, na.rm = TRUE)
           ))
  )

# View the descriptive statistics
print(descriptive_stats)

# Save the descriptive statistics as a CSV file
write.csv(descriptive_stats, file = "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/descriptive_stats.csv", row.names = FALSE)

install.packages("knitr")
install.packages("kableExtra")
library(knitr)
library(kableExtra)

# Filter data for men only and exclude rows with NA ages
men_data <- combined_data_complete_cases %>%
  filter(RIAGENDR == 1 & !is.na(RIDAGEYR))

# Define a function to calculate unweighted descriptive statistics for men
calculate_unweighted_stats_men <- function(data) {
  data %>%
    group_by(RIDRETH3, age_group) %>%
    summarise(
      n = n(),
      Mean_Age = mean(RIDAGEYR, na.rm = TRUE),
      SD_Age = sd(RIDAGEYR, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Generate unweighted descriptive statistics for men
unweighted_stats_men <- calculate_unweighted_stats_men(men_data)

# Decode race for readability
unweighted_stats_men <- unweighted_stats_men %>%
  mutate(
    race = recode(RIDRETH3,
                  `1` = "Mexican American",
                  `2` = "Other Hispanic",
                  `3` = "Non-Hispanic White",
                  `4` = "Non-Hispanic Black",
                  `6` = "Non-Hispanic Asian",
                  `7` = "Other Race - Including Multi-Racial")
  )

# Select and reorder columns for the table
unweighted_stats_men <- unweighted_stats_men %>%
  select(race, age_group, n, Mean_Age, SD_Age)

# Create and format the table
unweighted_stats_men %>%
  kable("html", col.names = c("Race", "Age Group", "N", "Mean Age", "SD Age"), 
        caption = "Unweighted Descriptive Statistics for Men by Race and Age Group") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  column_spec(1:5, width = "2em") %>%
  save_kable("/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/unweighted_descriptive_stats_men_table.html")

# Filter data for men only and exclude rows with NA ages
men_data <- combined_data_complete_cases %>%
  filter(RIAGENDR == 1 & !is.na(RIDAGEYR))

# Create a survey design object for men
nhanes_design_men <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC4YR,
  data = men_data,
  nest = TRUE
)

# Define age groups
men_data <- men_data %>%
  mutate(age_group = case_when(
    RIDAGEYR >= 20 & RIDAGEYR <= 39 ~ "20-39",
    RIDAGEYR >= 40 & RIDAGEYR <= 59 ~ "40-59",
    RIDAGEYR >= 60 ~ "60+",
    TRUE ~ NA_character_ 
  ))

# Calculate weighted n (number of participants)
weighted_n <- svyby(~one, ~RIDRETH3 + age_group, nhanes_design_men, svytotal, na.rm = TRUE)
weighted_n <- weighted_n %>%
  rename(n = one) %>%
  mutate(n = as.numeric(n))

# Decode race for readability
weighted_n <- weighted_n %>%
  mutate(
    race = recode(RIDRETH3,
                  `1` = "Mexican American",
                  `2` = "Other Hispanic",
                  `3` = "Non-Hispanic White",
                  `4` = "Non-Hispanic Black",
                  `6` = "Non-Hispanic Asian",
                  `7` = "Other Race - Including Multi-Racial")
  )

# Select and reorder columns for the table
weighted_n <- weighted_n %>%
  select(race, age_group, n)

# Print the table
print(weighted_n)

install.packages("webshot2")
library(webshot2)

# Create and format the table
weighted_n_table <- weighted_n %>%
  kable("html", col.names = c("Race", "Age Group", "Weighted N"), 
        caption = "Weighted Number of Participants for Men by Race and Age Group") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  column_spec(1:3, width = "3em")

# Save the table as a PDF
pdf_file_path <- "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/weighted_n_table.pdf"
save_kable(weighted_n_table, pdf_file_path)

