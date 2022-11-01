library(SummarizedExperiment)
library(sesame)
library(Rtsne)
library(stats)
library(factoextra)
library(RNAseqQC)
library(dplyr)

# Read in data as SummarizedExperiment
adni <- readRDS("~/Projects/havana/data/ADNI/adni.rds")

# Extract information
betas <- assays(adni)[[1]]
qcs <- metadata(adni)
adniQC <- do.call(rbind, lapply(qcs, as.data.frame))

# Density Plots
plot(density(adniQC$frac_dt))

# PCA
# Compute PCA
pca_res <- plot_pca(adni, n_feats = 500, na_frac = 0, show_plot = FALSE)

# Visualize by Variables of Interest
adni$PTGENDER <- as.factor(adni$PTGENDER)
plot_pca_scatters(adni,n_PCs = 4,color_by = "PTGENDER")

adni$PlateNumber <- as.factor(adni$PlateNumber)
plot_pca_scatters(adni,n_PCs = 4,color_by = "PlateNumber")

# Plot loadings
plot_loadings(pca_res, PC = 1)

# Perform TSNE 
res <- Rtsne(pca$x,pca=F,max_iter=2500,theta=0,verbose=T)

# Additional QCs
idat_dir_adni <- "/Users/sreepada/Library/CloudStorage/Box-Box/DataRepository/ADNI/Methylation/ADNI_iDAT_files/"
qcs_intensity <- openSesame(idat_dir_adni, prep="", func=sesameQC_calcStats, funs="intensity")

