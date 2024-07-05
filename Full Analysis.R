install.packages("dplyr")
install.packages("tidyr")
install.packages("survey")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("webshot2")
install.packages("purrr")
install.packages("knitr")
install.packages("kableExtra")
library(dplyr)
library(tidyr)
library(survey)
library(ggplot2)
library(gridExtra)
library(webshot2)
library(purrr)
library(knitr)
library(kableExtra)

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

#Waist over 35 in Women and 40 in Men
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BMX_H.XPT", tf <- tempfile(), mode="wb")
BMX_H <- foreign::read.xport(tf)[,c("SEQN", "BMIWAIST")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BMX_I.XPT", tf <- tempfile(), mode="wb")
BMX_I <- foreign::read.xport(tf)[,c("SEQN", "BMIWAIST")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/P_BMX.XPT", tf <- tempfile(), mode="wb")
BMX_P <- foreign::read.xport(tf)[,c("SEQN", "BMIWAIST")]

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

# Load necessary libraries
library(dplyr)
library(survey)

# Define age groups
combined_data <- combined_data %>%
  filter(!is.na(RIDAGEYR)) %>%  # Remove rows with NA ages
  mutate(Age_group = case_when(
    RIDAGEYR >= 0 & RIDAGEYR <= 19 ~ "0-19",
    RIDAGEYR >= 20 & RIDAGEYR <= 39 ~ "20-39",
    RIDAGEYR >= 40 & RIDAGEYR <= 59 ~ "40-59",
    RIDAGEYR >= 60 ~ "60+",
    TRUE ~ NA_character_
  ))

# Remove rows with NA Age_group
combined_data <- combined_data %>%
  filter(!is.na(Age_group))

# Group by race, gender, and age group
grouped_data <- combined_data %>%
  group_by(RIDRETH3, RIAGENDR, Age_group)

# View the first few rows of the grouped data
head(grouped_data)

# Create a summary table with all CBC variables
summary_table <- grouped_data %>%
  summarise(
    Count = n(),
    Median_LBXWBCSI = median(LBXWBCSI, na.rm = TRUE),
    Median_LBXLYPCT = median(LBXLYPCT, na.rm = TRUE),
    Median_LBXMOPCT = median(LBXMOPCT, na.rm = TRUE),
    Median_LBXNEPCT = median(LBXNEPCT, na.rm = TRUE),
    Median_LBXEOPCT = median(LBXEOPCT, na.rm = TRUE),
    Median_LBXBAPCT = median(LBXBAPCT, na.rm = TRUE),
    Median_LBXRBCSI = median(LBXRBCSI, na.rm = TRUE),
    Median_LBXHGB = median(LBXHGB, na.rm = TRUE),
    Median_LBXHCT = median(LBXHCT, na.rm = TRUE),
    Median_LBXMCVSI = median(LBXMCVSI, na.rm = TRUE),
    Median_LBXMCHSI = median(LBXMCHSI, na.rm = TRUE),
    Median_LBXMC = median(LBXMC, na.rm = TRUE),
    Median_LBXRDW = median(LBXRDW, na.rm = TRUE),
    Median_LBXPLTSI = median(LBXPLTSI, na.rm = TRUE),
    Median_LBXMPSI = median(LBXMPSI, na.rm = TRUE)
  )

# View the summary table
print(summary_table)

# Optional: create a survey design object with the grouped data
nhanes_design_grouped <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC4YR,
  data = grouped_data,
  nest = TRUE
)

# Check the summary of the survey design object for the grouped data
summary(nhanes_design_grouped)

# Check the structure and summary of the combined dataset
str(combined_data)
summary(combined_data)

# Check for missing values in important columns
sapply(combined_data, function(x) sum(is.na(x)))


# Function to print number of rows and sample data after each exclusion step
print_status <- function(data, step) {
  cat("\nNumber of rows after", step, ":", nrow(data), "\n")
  head(data)
}

# Check initial dimensions and summary
print("Initial data dimensions and summary:")
dim(combined_data)
summary(combined_data)

# Apply exclusion criteria: A1C >= 6.5%
combined_data_excluded_a1c <- combined_data %>%
  filter(LBXGH < 6.5 | is.na(LBXGH))
print_status(combined_data_excluded_a1c, "excluding A1C >= 6.5%")

# Create a new survey design object using the filtered data
nhanes_design_excluded_a1c <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_a1c, nest = TRUE)

# Apply exclusion criteria: Serum creatinine > 2.5 mg/dL
combined_data_excluded_creatinine <- combined_data_excluded_a1c %>%
  filter(LBXSCR < 2.5 | is.na(LBXSCR))
print_status(combined_data_excluded_creatinine, "excluding Serum creatinine > 2.5 mg/dL")

# Create a new survey design object using the filtered data
nhanes_design_excluded_creatinine <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_creatinine, nest = TRUE)

# Apply exclusion criteria: Hepatitis C positivity (LBXHCR == 1)
combined_data_excluded_hepc <- combined_data_excluded_creatinine %>%
  filter(LBXHCR != 1 | is.na(LBXHCR))
print_status(combined_data_excluded_hepc, "excluding Hepatitis C positivity")

# Create a new survey design object using the filtered data
nhanes_design_excluded_hepc <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_hepc, nest = TRUE)

# Apply exclusion criteria: Hepatitis B surface antigen (LBDHBG == 1) or core antibody (LBXHBC == 1)
combined_data_excluded_hepb <- combined_data_excluded_hepc %>%
  filter((LBDHBG != 1 | is.na(LBDHBG)) & (LBXHBC != 1 | is.na(LBXHBC)))
print_status(combined_data_excluded_hepb, "excluding Hepatitis B surface antigen or core antibody")

# Create a new survey design object using the filtered data
nhanes_design_excluded_hepb <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_hepb, nest = TRUE)

# Apply exclusion criteria for 2015-2016 data: ALQ141Q >= 104
combined_data_excluded_alcohol_2015_2016 <- combined_data_excluded_hepb %>%
  filter(is.na(ALQ141Q) | ALQ141Q < 104)
print_status(combined_data_excluded_alcohol_2015_2016, "excluding ALQ141Q >= 104 (2015-2016)")

# Apply exclusion criteria for 2017-2018 data: ALQ142 in 1, 2, 3, or 4
combined_data_excluded_alcohol <- combined_data_excluded_alcohol_2015_2016 %>%
  filter(!(ALQ142 %in% c(1, 2, 3, 4)))
print_status(combined_data_excluded_alcohol, "excluding ALQ142 in 1, 2, 3, or 4 (2017-2018)")

# Create a new survey design object using the filtered data
nhanes_design_excluded_alcohol <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_alcohol, nest = TRUE)

# Apply exclusion criteria: Pregnancy (RHD143 == 1) and Breastfeeding (RHQ200 == 1) for females
combined_data_excluded_pregnancy_breastfeeding <- combined_data_excluded_alcohol %>%
  filter(!(RIAGENDR == 2 & (RHD143 == 1 | RHQ200 == 1)))
print_status(combined_data_excluded_pregnancy_breastfeeding, "excluding Pregnancy and Breastfeeding")

# Create a new survey design object using the filtered data
nhanes_design_excluded_pregnancy_breastfeeding <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_pregnancy_breastfeeding, nest = TRUE)

# Apply exclusion criteria: BMI (BMIWAIST >35 in Women and >40 in Men)
combined_data_excluded_bmi <- combined_data_excluded_pregnancy_breastfeeding %>%
  filter((RIAGENDR == 2 & (BMIWAIST <= 35 | is.na(BMIWAIST))) | 
           (RIAGENDR == 1 & (BMIWAIST <= 40 | is.na(BMIWAIST))))
print_status(combined_data_excluded_bmi, "excluding BMI (BMIWAIST >35 in Women and >40 in Men)")

# Create a new survey design object using filtered data
nhanes_design_excluded_bmi <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_bmi, nest = TRUE)

# Apply exclusion criteria: Age < 20
combined_data_excluded_age <- combined_data_excluded_bmi %>%
  filter(RIDAGEYR > 20)
print_status(combined_data_excluded_age, "excluding Age < 20")

# Create a new survey design object using the filtered data
nhanes_design_excluded_age <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_excluded_age, nest = TRUE)

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
print_status(combined_data_complete_cases, "excluding Missing CBCs")

# Ensure age groups are defined correctly after exclusions
combined_data_complete_cases <- combined_data_complete_cases %>%
  filter(!is.na(RIDAGEYR)) %>%  # Remove rows with NA ages
  mutate(Age_group = case_when(
    RIDAGEYR >= 0 & RIDAGEYR <= 19 ~ "0-19",
    RIDAGEYR >= 20 & RIDAGEYR <= 39 ~ "20-39",
    RIDAGEYR >= 40 & RIDAGEYR <= 59 ~ "40-59",
    RIDAGEYR >= 60 ~ "60+",
    TRUE ~ NA_character_
  ))

# Remove rows with NA Age_group
combined_data_complete_cases <- combined_data_complete_cases %>%
  filter(!is.na(Age_group))

# Create a survey design object using the filtered data
nhanes_design_complete_cases <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR, data = combined_data_complete_cases, nest = TRUE)

# Function to filter data by race, gender, and age group
filter_data <- function(data, race, gender, age_group) {
  data %>%
    filter(RIDRETH3 == race & RIAGENDR == gender & Age_group == age_group)
}

races <- c(1, 2, 3, 4, 6)  # Exclude race code 7 (other, mixed race)
genders <- c(1, 2)  # 1 for Male, 2 for Female
age_groups <- c("20-39", "40-59", "60+")

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

# Remove cases where age_group is NA and race is "other, mixed race"
combined_data_complete_cases <- combined_data_complete_cases %>%
  filter(!is.na(Age_group) & RIDRETH3 != 7)

# Count the number of cases for each combination of gender, race, and age group
count_cases <- combined_data_complete_cases %>%
  group_by(RIDRETH3, RIAGENDR, Age_group) %>%
  summarise(count = n()) %>%
  ungroup()

# Print the count cases
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
  group_by(RIDRETH3, RIAGENDR, Age_group) %>%
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

# Filter data for men only and exclude rows with NA ages
men_data <- combined_data_complete_cases %>%
  filter(RIAGENDR == 1 & !is.na(RIDAGEYR))

# Define a function to calculate unweighted descriptive statistics for men
calculate_unweighted_stats_men <- function(data) {
  data %>%
    group_by(RIDRETH3, Age_group) %>%
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
  select(race, Age_group, n, Mean_Age, SD_Age)

# Create and format the table
unweighted_stats_men %>%
  kable("html", col.names = c("Race", "Age Group", "N", "Mean Age", "SD Age"), 
        caption = "Unweighted Descriptive Statistics for Men by Race and Age Group") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  column_spec(1:5, width = "2em") %>%
  save_kable("/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/unweighted_descriptive_stats_men_table.html")


# Install and load necessary packages
install.packages("webshot")
library(webshot)

# Install PhantomJS for webshot2
webshot::install_phantomjs()

# Filter data for men only, exclude rows with NA ages and "other race" category
men_data <- combined_data_complete_cases %>%
  filter(RIAGENDR == 1 & !is.na(RIDAGEYR) & RIDRETH3 != 7)

# Add a constant variable `one` with value 1
men_data <- men_data %>%
  mutate(one = 1)

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
  mutate(Age_group = case_when(
    RIDAGEYR >= 20 & RIDAGEYR <= 39 ~ "20-39",
    RIDAGEYR >= 40 & RIDAGEYR <= 59 ~ "40-59",
    RIDAGEYR >= 60 ~ "60+",
    TRUE ~ NA_character_ 
  ))

# Calculate weighted n (number of participants)
weighted_n <- svyby(~one, ~RIDRETH3 + Age_group, nhanes_design_men, svytotal, na.rm = TRUE)
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
                  `6` = "Non-Hispanic Asian")
  )

# Select and reorder columns for the table
weighted_n <- weighted_n %>%
  select(race, Age_group, n)

# Print the table
print(weighted_n)

# Create and format the table
weighted_n_table <- weighted_n %>%
  kable("html", col.names = c("Race", "Age Group", "Weighted N"), 
        caption = "Weighted Number of Participants for Men by Race and Age Group") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  column_spec(1:3, width = "3em")

# Save the table as a PDF
pdf_file_path <- "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/weighted_n_table.pdf"
save_kable(weighted_n_table, pdf_file_path)

#Calculating non-parametric reference ranges for the CBC parameters
## Define a function to calculate weighted percentiles
calculate_weighted_percentiles <- function(variable, design) {
  result <- svyby(
    ~ get(variable),
    ~ RIDRETH3 + Age_group,
    design,
    svyquantile,
    quantiles = c(0.025, 0.5, 0.975),
    ci = TRUE,
    na.rm = TRUE
  )
  result$Variable <- variable
  return(result)
}

# List of CBC variables
cbc_vars <- c("LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", 
              "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", 
              "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")

# Calculate percentiles for each CBC variable
percentile_list <- lapply(cbc_vars, calculate_weighted_percentiles, nhanes_design_men)

# Combine the results into a single data frame
percentile_df <- do.call(rbind, percentile_list)

# Rename the columns for readability
colnames(percentile_df) <- c("Race", "Age Group", "2.5th Percentile", "50th Percentile", "97.5th Percentile", "Variable")

# Decode race for readability
percentile_df$Race <- recode(percentile_df$Race,
                             `1` = "Mexican American",
                             `2` = "Other Hispanic",
                             `3` = "Non-Hispanic White",
                             `4` = "Non-Hispanic Black",
                             `6` = "Non-Hispanic Asian")

# Print the percentiles
print(percentile_df)

# Save the percentiles as a CSV file
write.csv(percentile_df, file = "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/weighted_percentiles_with_variables.csv", row.names = FALSE)


#KW Tests

# Load necessary packages
library(dplyr)
library(stats)
library(FSA)
library(ggplot2)

# Function to run Kruskal-Wallis test based on filtered data for men
run_kruskal_wallis_actual_data <- function(data, variable, race) {
  # Filter data for the specified CBC variable, race, and excluding age group 0-19
  test_data <- data %>%
    filter(RIDRETH3 == race & Age_group != "0-19" & !is.na(.data[[variable]])) %>%
    mutate(Age_group = as.factor(Age_group)) # Ensure Age_group is a factor
  
  # Print the first few rows of the filtered data to ensure it's correct
  print(paste("Filtered data for variable", variable, "and race", race))
  print(head(test_data))
  
  # Check if there are at least two groups to compare
  if (n_distinct(test_data$Age_group) > 1) {
    # Perform the Kruskal-Wallis test
    kw_test <- kruskal.test(as.formula(paste(variable, "~ Age_group")), data = test_data)
    
    # Print the results
    print(paste("Kruskal-Wallis test for", variable, "in race", race))
    print(kw_test)
    
    # Return the result as a data frame and the filtered data
    result <- data.frame(
      Variable = variable,
      Race = race,
      p_value = kw_test$p.value
    )
    
    return(list(result = result, test_data = test_data))
  } else {
    print("Not enough groups to perform the Kruskal-Wallis test.")
    return(NULL)
  }
}

# Manually enter the CBC variable and race
variable <- "LBXPLTSI" # Replace with the desired CBC variable
race <- 4 # Replace with the desired race code (e.g., 6 for Non-Hispanic Asian)

# Run the Kruskal-Wallis test using the men_data dataset
kw_result <- run_kruskal_wallis_actual_data(men_data, variable, race)

# Print the result if not null
if (!is.null(kw_result)) {
  print(kw_result$result)
  
  # Extract the test data for further analysis
  test_data <- kw_result$test_data
  
  # Perform Dunn's test for post-hoc analysis
  dunn_result <- dunnTest(as.formula(paste(variable, "~ Age_group")), data = test_data, method = "bonferroni")
  
  # Print the Dunn's test result
  print(dunn_result)
  
  # Visualize the data with a boxplot
  ggplot(test_data, aes(x = Age_group, y = .data[[variable]], fill = Age_group)) +
    geom_boxplot() +
    labs(title = paste("Distribution of", variable, "by Age Group within Non-Hispanic Asian Race (Males only)"),
         x = "Age Group",
         y = variable) +
    theme_minimal()
} else {
  print("Kruskal-Wallis test result is NULL.")
}

# Across Race by Age
# Load necessary packages
library(dplyr)
library(stats)
library(FSA)
library(ggplot2)

# Function to run Kruskal-Wallis test and Dunn's test for each age group across races
run_kruskal_wallis_by_age_group <- function(data, variable) {
  age_groups <- c("20-39", "40-59", "60+")
  
  results <- list()
  
  for (age_group in age_groups) {
    # Filter data for the specified age group and CBC variable
    test_data <- data %>%
      filter(Age_group == age_group & !is.na(.data[[variable]])) %>%
      mutate(Race = factor(RIDRETH3, labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian")))
    
    # Print the first few rows of the filtered data to ensure it's correct
    print(paste("Filtered data for variable", variable, "and age group", age_group))
    print(head(test_data))
    
    # Check if there are at least two groups to compare
    if (n_distinct(test_data$Race) > 1) {
      # Perform the Kruskal-Wallis test
      kw_test <- kruskal.test(as.formula(paste(variable, "~ Race")), data = test_data)
      
      # Print the results
      print(paste("Kruskal-Wallis test for", variable, "in age group", age_group))
      print(kw_test)
      
      # Store the result
      results[[age_group]] <- list(kw_test = kw_test)
      
      # If significant, perform Dunn's test
      if (kw_test$p.value < 0.05) {
        dunn_result <- dunnTest(as.formula(paste(variable, "~ Race")), data = test_data, method = "bonferroni")
        
        # Print the Dunn's test result
        print(dunn_result)
        
        # Store the Dunn's test result
        results[[age_group]]$dunn_result <- dunn_result
        
        # Visualize the data with a boxplot
        p <- ggplot(test_data, aes(x = Race, y = .data[[variable]], fill = Race)) +
          geom_boxplot() +
          labs(title = paste("Distribution of", variable, "by Race within Age Group", age_group),
               x = "Race",
               y = variable) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
      }
    } else {
      print(paste("Not enough groups to perform the Kruskal-Wallis test for age group", age_group))
    }
  }
  
  return(results)
}

# List of CBC variables to analyze
cbc_vars <- c("LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT", "LBXEOPCT", 
              "LBXBAPCT", "LBXRBCSI", "LBXHGB", "LBXHCT", "LBXMCVSI", 
              "LBXMCHSI", "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI")

# Define the path for the PDF file
pdf_path <- "/Users/seanchickery/Library/Mobile Documents/com~apple~CloudDocs/R NHANES CBC 2015 thru 2018/kruskal_wallis_results.pdf"

# Open PDF device
pdf(file = pdf_path)

# Run the Kruskal-Wallis test and Dunn's test for each CBC variable across age groups and save to PDF
results <- lapply(cbc_vars, run_kruskal_wallis_by_age_group, data = men_data)

# Close PDF device
dev.off()

# Optionally, you can print or save the results
print(results)

# Load necessary package
library(dplyr)

# Count the number of unique SEQNs in men_data
num_seqns <- men_data %>% 
  summarise(unique_seqns = n_distinct(SEQN))

# Print the result
print(paste("Number of unique SEQNs in men_data:", num_seqns$unique_seqns))



