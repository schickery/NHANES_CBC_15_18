install.packages("survey")
library(survey)
library(dplyr)
install.packages("ggplot2")
library(ggplot2)

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
HEPC_J <- foreign::read.xport(tf)[,c("SEQN", "LBDHCI")]
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

# Combine all datasets into one dataframe (assuming SEQN is the unique identifier across all datasets)
combined_data <- reduce(list(DEMO, CBC, BIOPRO, HEPBD, HEPC, HSQ, ALQ, SMQ, GHB, RHQ), full_join, by = "SEQN")

# Count the number of unique SEQN values
unique_seqn_count <- n_distinct(combined_data$SEQN)
print(paste("Number of unique SEQN in the combined dataset:", unique_seqn_count))

# Decoding gender and race if necessary
DEMO <- DEMO %>%
  mutate(
    gender = recode(RIAGENDR, `1` = "Male", `2` = "Female"),
    race = recode(RIDRETH3,
                  `1` = "Mexican American",
                  `2` = "Other Hispanic",
                  `3` = "Non-Hispanic White",
                  `4` = "Non-Hispanic Black",
                  `6` = "Non-Hispanic Asian",
                  `7` = "Other Race - Including Multi-Racial")
  )

# Summarizing counts by gender and race
gender_race_counts <- DEMO %>%
  group_by(gender, race) %>%
  summarise(count = n(), .groups = 'drop')

print(gender_race_counts)

#Create a Table for m/F by Race-----------------------------------
install.packages("tidyr")
install.packages("gridExtra")
install.packages("grid")
library(tidyr)
library(gridExtra)
library(grid)

# Summarizing counts by gender and race
gender_race_counts <- DEMO %>%
  group_by(gender, race) %>%
  summarise(count = n(), .groups = 'drop')

# Pivoting the data to wide format
gender_race_table <- gender_race_counts %>%
  pivot_wider(names_from = race, values_from = count, values_fill = list(count = 0))

print(gender_race_table)

# Convert the table to a format suitable for ggplot2
gender_race_table_long <- gender_race_table %>%
  pivot_longer(-gender, names_to = "race", values_to = "count")

# Create a table plot using ggplot2
table_plot <- ggplot(gender_race_table_long, aes(x = race, y = gender)) +
  geom_tile(aes(fill = count), color = "white") +
  geom_text(aes(label = count), color = "black") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Number of Males and Females by Race",
       fill = "Count",
       x = "Race",
       y = "Gender")

# Save the table plot as a PDF
ggsave("gender_race_table.pdf", plot = table_plot, width = 10, height = 6)
