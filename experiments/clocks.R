# Libraries
library(dplyr)
library(methylclock)
library(sesame)
library(SummarizedExperiment)
library(stats)
library(smplot)
library(ggplot2)
library(lubridate)

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
betasHorvath <- tibble::rownames_to_column(as.data.frame(betas),"ProbeID") # Get rownames (CpGs) as first column, Horvath format
betasHorvath <- subset(betasHorvath,ProbeID %in% coefHorvath$CpGmarker) # subset betas matrix to only include Horvath CpGs
temp <- DNAmAge(betasHorvath,clocks=c("Horvath")) 
adni$AgeHorvath <- temp$Horvath
  
# Compute Shireby CorticalAge
clockDir <- "CorticalClock/PredCorticalAge/"
source(paste0(clockDir,"RunCorticalClock.R"))

coefCortical <- read.table(paste0(clockDir,"CorticalClockCoefs.txt",sep=""),stringsAsFactor=F,header=T)
betasCortical <- as.data.frame(betas[rownames(betas) %in% coefCortical$probe,])
RunCorticalClock(betasCortical, as.data.frame(colData(adni)), clockDir, "barcodes", "Age",resultsdir = 'data/ADNI/',outfile = "ADNI")

corticalAge <- read.csv("data/ADNI/ADNI_CorticalPred.csv")
adni$AgeCortical <- corticalAge$brainpred

# Write coldata to csv for easy manipulation
write.csv(colData(adni),'data/ADNI/ClocksAndImaging.csv')
