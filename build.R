#!/usr/bin/env Rscript

install.packages("tidyr")
install.packages("devtools") 
library(devtools)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CytoML", ask=FALSE, force=TRUE)
BiocManager::install("cytolib", ask=FALSE, force=TRUE)

install_github("RGLab/flowWorkspace", ref="bioc_3.13", dependencies=TRUE, upgrade_dependencies=TRUE, ask = F)
