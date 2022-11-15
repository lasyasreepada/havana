library(ggplot2)
library(smplot)
library(dplyr)

# Read in ADNI Dataset

# Plotting
gapMedian <- summarise(group_by(as.data.frame(colData(adni)), as.factor(DX1)), medianHorvathGap = median(AgeHorvathGap))

# Correlation between Horvath and Age, stratified by DX
ggplot(as.data.frame(colData(adni)), aes(x=Age, y=AgeHorvath, color=as.factor(DX1))) +
  geom_point() +
  geom_smooth(aes(color=as.factor(DX1)),method="lm") +
  geom_abline(slope = 1, intercept = 0) + 
  sm_corr_theme() +
  sm_statCorr(line_type = 'blank')

ggplot(as.data.frame(colData(adni)), aes(x=Age, y=cortToMTL, color=as.factor(DX1))) +
  geom_point() +
  geom_smooth(aes(color=as.factor(DX1)),method="lm") +
  geom_abline(slope = 1, intercept = 0) + 
  sm_corr_theme() +
  sm_statCorr(line_type = 'blank')

# Boxplots of Age difference stratified by DX
ggplot(as.data.frame(colData(adni)), aes(x=as.factor(DX1), y=AgeHorvathGap, color=as.factor(DX1))) +
  geom_boxplot() 

# Boxplot by APOE
ggplot(as.data.frame(colData(adni)), aes(x=as.factor(APOE), y=AgeHorvathGap, color=as.factor(APOE))) +
  geom_boxplot()

# Boxplot by APOE
ggplot(as.data.frame(colData(adni)), aes(x=as.factor(DX1), y=cortToMTL, color=as.factor(DX1))) +
  geom_boxplot()