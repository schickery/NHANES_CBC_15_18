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
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT", tf <- tempfile(), mode="wb")
DEMO_I <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDEXPRG", "RIDRETH3", "SDMVSTRA","SDMVPSU","WTMEC2YR")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT", tf <- tempfile(), mode="wb")
DEMO_J <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDEXPRG", "RIDRETH3", "SDMVSTRA","SDMVPSU","WTMEC2YR")]


#CBC Data (CBC)
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/CBC_I.XPT", tf <- tempfile(), mode="wb")
CBC_I <- foreign::read.xport(tf)[,c("SEQN", "LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/CBC_J.XPT", tf <- tempfile(), mode="wb")
CBC_J <- foreign::read.xport(tf)[,c("SEQN", "LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")]

#Biochem for Creatinine and Glucose
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BIOPRO_I.XPT", tf <- tempfile(), mode="wb")
BIOPRO_I <- foreign::read.xport(tf)[,c("SEQN", "LBXSCR", "LBXSGL")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BIOPRO_J.XPT", tf <- tempfile(), mode="wb")
BIOPRO_J <- foreign::read.xport(tf)[,c("SEQN", "LBXSCR", "LBXSGL")]

#Hepatitis B Surface Ag and Core Ab
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HEPBD_I.XPT", tf <- tempfile(), mode="wb")
HEPBD_I <- foreign::read.xport(tf)[,c("SEQN", "LBXHBC", "LBDHBG")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPBD_J.XPT", tf <- tempfile(), mode="wb")
HEPBD_J <- foreign::read.xport(tf)[,c("SEQN", "LBXHBC", "LBDHBG")]

#HEP C Confirmed Data
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPC_J.XPT", tf <- tempfile(), mode="wb")
HEPC_J <- foreign::read.xport(tf)[,c("SEQN", "LBXHCR")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HEPC_I.XPT", tf <- tempfile(), mode="wb")
HEPC_I <- foreign::read.xport(tf)[,c("SEQN", "LBXHCR")]

#Blood Donation & Illness (Exclude Past 12 weeks donators)
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSQ_I.XPT", tf <- tempfile(), mode="wb")
HSQ_I <- foreign::read.xport(tf)[,c("SEQN", "HSQ500", "HSQ510", "HSQ520","HSQ571", "HSQ580")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HSQ_J.XPT", tf <- tempfile(), mode="wb")
HSQ_J <- foreign::read.xport(tf)[,c("SEQN", "HSQ500", "HSQ510", "HSQ520","HSQ571", "HSQ580")]

#Alcohol Use
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/ALQ_I.XPT", tf <- tempfile(), mode="wb")
ALQ_I <- foreign::read.xport(tf)[,c("SEQN", "ALQ141Q", "ALQ141U")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/ALQ_J.XPT", tf <- tempfile(), mode="wb")
ALQ_J <- foreign::read.xport(tf)[,c("SEQN", "ALQ142")]

#Smoking
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SMQ_I.XPT", tf <- tempfile(), mode="wb")
SMQ_I <- foreign::read.xport(tf)[,c("SEQN", "SMQ040")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SMQ_J.XPT", tf <- tempfile(), mode="wb")
SMQ_J <- foreign::read.xport(tf)[,c("SEQN", "SMQ040")]

#A1C
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/GHB_I.XPT", tf <- tempfile(), mode="wb")
GHB_I <- foreign::read.xport(tf)[,c("SEQN", "LBXGH")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/GHB_J.XPT", tf <- tempfile(), mode="wb")
GHB_J <- foreign::read.xport(tf)[,c("SEQN", "LBXGH")]

#Pregnancy and Breast Feeding
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/RHQ_I.XPT", tf <- tempfile(), mode="wb")
RHQ_I <- foreign::read.xport(tf)[,c("SEQN", "RHD143", "RHQ200")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/RHQ_J.XPT", tf <- tempfile(), mode="wb")
RHQ_J <- foreign::read.xport(tf)[,c("SEQN", "RHD143", "RHQ200")]

#Append Files
DEMO <- bind_rows(DEMO_I, DEMO_J)
CBC <- bind_rows(CBC_I, CBC_J)
BIOPRO <- bind_rows(BIOPRO_I, BIOPRO_J)
HEPBD <- bind_rows(HEPBD_I, HEPBD_J)
HEPC <- bind_rows(HEPC_I, HEPC_J)
HSQ <- bind_rows(HSQ_I, HSQ_J)
ALQ <- bind_rows(ALQ_I, ALQ_J)
SMQ <- bind_rows(SMQ_I, SMQ_J)
GHB <- bind_rows(GHB_I, GHB_J)
RHQ <- bind_rows(RHQ_I, RHQ_J)

#------------------------------------------
#Determine number of SEQNs before exclusions

install.packages("purrr")
library(purrr)

# Adjust weights for 4 years of data
DEMO <- DEMO %>%
  mutate(WTMEC4YR = WTMEC2YR / 2)

# Combine all datasets into one dataframe
combined_data <- reduce(list(DEMO, CBC, BIOPRO, HEPBD, HEPC, HSQ, ALQ, SMQ, GHB, RHQ), full_join, by = "SEQN")

# Check that the necessary columns exist
required_columns <- c("SEQN", "SDMVPSU", "SDMVSTRA", "WTMEC4YR")
missing_columns <- setdiff(required_columns, colnames(combined_data))

if(length(missing_columns) > 0) {
  stop(paste("Missing columns:", paste(missing_columns, collapse = ", ")))
}

# Create a survey design object using the combined data
nhanes_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data, nest = TRUE)


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
combined_data_excluded_blood_donation <- combined_data_excluded_hepb %>%
  filter(!(HSQ580 %in% c(1, 2, 3)))

# Create a new survey design object using the filtered data
nhanes_design_excluded_blood_donation <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_blood_donation, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_blood_donation, nhanes_design_excluded_blood_donation, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_blood_donation.pdf")

# Apply exclusion criteria for 2015-2016 data: ALQ141Q >= 54----------------------------------------------
combined_data_excluded_alcohol_2015_2016 <- combined_data_excluded_blood_donation %>%
  filter((is.na(ALQ141Q) | ALQ141Q < 54))

# Apply exclusion criteria for 2017-2018 data: ALQ142 in 1, 2, 3, or 4
combined_data_excluded_alcohol <- combined_data_excluded_alcohol_2015_2016 %>%
  filter(!(ALQ142 %in% c(1, 2, 3, 4)))

# Create a new survey design object using the filtered data
nhanes_design_excluded_alcohol <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_alcohol, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_alcohol, nhanes_design_excluded_alcohol, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_alcohol.pdf")

# Apply exclusion criteria: Pregnancy (RHD143 == 1) and Breastfeeding (RHQ200 == 1) for females
combined_data_excluded_pregnancy_breastfeeding <- combined_data_excluded_alcohol %>%
  filter(!(RIAGENDR == 2 & (RHD143 == 1 | RHQ200 == 1)))

# Create a new survey design object using the filtered data
nhanes_design_excluded_pregnancy_breastfeeding <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_pregnancy_breastfeeding, nest = TRUE)

# Summarize counts after exclusion and save as PDF
summarize_counts(combined_data_excluded_pregnancy_breastfeeding, nhanes_design_excluded_pregnancy_breastfeeding, "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/summary_after_exclusion_pregnancy_breastfeeding.pdf")





