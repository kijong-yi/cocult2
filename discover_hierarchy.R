library(tidyverse)
library(iNEXT)
library(vegan)
library(doMC)
library(Hmisc)
library(ComplexHeatmap)
registerDoMC(24)


source("~/Projects/cocult/cocult2/Wp.R")

A.std2 = kjyi::read_rdata("~/Projects/cocult/cocult2/snapshot_20191008.RData", "A.std2")
Data = t(A.std2)[,sample(66579,1000)]
dim(Data)
res = Wp_R(Data, p = 1.5, RequiredK = 1L, trace = F)


