library(dplyr)
library(sesame)
library(snow)
library(stats)
library(RSQLite)
library(BiocParallel)
library(SummarizedExperiment)

# DATA IO
idat_dir_adni <- "/Users/sreepada/Library/CloudStorage/Box-Box/DataRepository/ADNI/Methylation/ADNI_iDAT_files/"
idat_dir_abc <- "/Users/sreepada/Library/CloudStorage/Box-Box/DataRepository/ABC/Methylation/"

# ADNI
adni_betas <- openSesame(idat_dir_adni, func = getBetas)
adni_qcs <- openSesame(idat_dir_adni, prep="", func=sesameQC_calcStats, funs="detection")

# ABC
abc_betas <- openSesame(idat_dir_abc, func = getBetas)
abc_qcs <- openSesame(idat_dir_abc, prep="", func=sesameQC_calcStats, funs="detection")

# Get names
adni_samples <- colnames(adni_betas)
abc_samples <- colnames(abc_betas)

# ABC
# Read in subject information (coldata)
abc_m_file <- "data/ABC/ABC_Methylation_Samples.csv"
abc_m_samples <- read.csv(abc_m_file) # ABC

# Create new column 'barcodes' for matching methylation data samples with beta matrix colnames in ABC
abc_m_samples$barcodes <- paste(abc_m_samples$Sentrix_ID, "_", abc_m_samples$Sentrix_Position)
abc_m_samples$barcodes <- str_replace_all(abc_m_samples$barcodes, pattern=" ", repl="")

# Match methylation samples with annotation based on 'barcodes'
abc_colData <- abc_m_samples[match(abc_samples,abc_m_samples$barcodes),]
rownames(abc_colData) <- abc_colData$barcodes

# Save as SummarizedExperiment
abc_se <- SummarizedExperiment(assays = abc_betas, colData=abc_colData, metadata = abc_qcs)

# Save to RDS file
saveRDS(abc_se, "data/ABC/methylation.rds")

# ADNI
# Read in subject information (coldata)
adni_m_file <- "data/ADNI/ADNI_DNA_Methylation_SampleAnnotation_20170530.xlsx"
adni_m_samples <- readxl::read_xlsx(adni_m_file) # ADNI
adni_samples <- as.data.frame(adni_samples)
colnames(adni_samples) <- "barcodes"

# Match methylation samples with annotation based on 'barcodes'
adni_colData <- as.data.frame(full_join(adni_samples, adni_m_samples, by="barcodes"))

# Save as SummarizedExperiment
adni_se <- SummarizedExperiment(assays = adni_betas, colData=adni_colData, metadata = adni_qcs)

# Save to RDS file
saveRDS(adni_se, "data/ADNI/methylation.rds")

# Some random stuff
# https://www.bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/inferences.html
# identical(colnames(assay(se)), rownames(colData(se))

