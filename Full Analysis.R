install.packages("survey")
library(survey)
library(dplyr)
install.packages("ggplot2")
library(ggplot2)

#Install Files from NHANES Database
#Demographics files
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT", tf <- tempfile(), mode="wb")
DEMO_J <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDEXPRG", "RIDRETH3", "SDMVSTRA","SDMVPSU","WTMEC2YR")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT", tf <- tempfile(), mode="wb")
DEMO_I <- foreign::read.xport(tf)[,c("SEQN","RIAGENDR","RIDAGEYR","RIDEXPRG", "RIDRETH3", "SDMVSTRA","SDMVPSU","WTMEC2YR")]

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
HEPC_J <- foreign::read.xport(tf)[,c("SEQN", "LBDHCI")]

#Blood Donation (Exclude Past 12 weeks donators)
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSQ_I.XPT", tf <- tempfile(), mode="wb")
HSQ_I <- foreign::read.xport(tf)[,c("SEQN", "HSQ571", "HSQ580")]
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HSQ_J.XPT", tf <- tempfile(), mode="wb")
HSQ_J <- foreign::read.xport(tf)[,c("SEQN", "HSQ571", "HSQ580")]

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
