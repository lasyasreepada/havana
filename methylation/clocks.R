# Libraries
library(dplyr)
library(methylclock)
library(sesame)
library(SummarizedExperiment)
library(stats)

# DATA IO
# Read in data as SummarizedExperiment
adni <- readRDS("~/Projects/havana/data/ADNI/adni.rds")

# Extract information
betas <- assays(adni)[[1]]

# Check missing CpGs by clock
cpgs.missing <- checkClocks(betas)

# Compute Horvath
betas2 <- tibble::rownames_to_column(as.data.frame(betas),"ProbeID") # Get rownames (CpGs) as first column, Horvath format
adni$HorvathAge <- DNAmAge(betas2,clocks=c("Horvath"), cell.count = FALSE) 