# Libraries
library(dplyr)
library(methylclock)
library(sesame)
library(SummarizedExperiment)
library(stats)
library(smplot)

# DATA IO
# Read in data as SummarizedExperiment
adni <- readRDS("~/Projects/havana/data/ADNI/adni.rds")

# Compute Chronological Age
adni$Age <- time_length(difftime(as.Date(adni$SCANDATE, format = "%Y-%m-%d"), as.Date(adni$birthdate, format = "%m/%d/%Y")),"years")

# Extract information
betas <- assays(adni)[[1]]

# Check missing CpGs by clock
cpgs.missing <- checkClocks(betas)

# Compute Horvath Age
betas2 <- tibble::rownames_to_column(as.data.frame(betas),"ProbeID") # Get rownames (CpGs) as first column, Horvath format
betasHorvath <- subset(betas2,ProbeID %in% coefHorvath$CpGmarker) # subset betas matrix to only include Horvath CpGs
temp <- DNAmAge(betasHorvath,clocks=c("Horvath")) 
adni$AgeHorvath <- temp$Horvath

# Compute Shireby CorticalAge
clockDir <- "CorticalClock/PredCorticalAge/"
source(path(clockDir, "CorticalClock.R"))
CorticalClock(betas, as.data.frame(colData(adni)), clockDir, "barcodes", "Age")

# Plotting
ggplot(as.data.frame(colData(adni)), aes(x=Age, y=AgeHorvath, color=as.factor(DX1))) +
  geom_point() +
  geom_smooth(aes(color=as.factor(DX1)),method="lm") +
  geom_abline(slope = 1, intercept = 0) + 
  sm_corr_theme()
