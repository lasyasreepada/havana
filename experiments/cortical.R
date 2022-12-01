library(dplyr)

adniCort <- read.csv("data/ADNI/ADNI_thickness_2021.csv")

cortical <- list(names(adniCort))
write.csv(cortical,"data/ADNI/corticalRegions.csv")