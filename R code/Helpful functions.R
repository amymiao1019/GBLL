# Need to library all the packages used:
library(ggplot2)
library(tidyr)
library(demography)
library(tseries)
library(gridExtra)
library(cowplot)
library(MortCast)
library(hwwntest)
library(tidyverse)
library(forecast)
library(lubridate)
library(lmtest)
library(dplyr)
library(MCS)
library(purrr)
library(cluster)
library(readxl)
library(patchwork)
library(reshape2)
library(pheatmap)
library(maps)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(wbstats)
library(ggrepel)
library(panelvar)
library(plm)
library(ggpubr)




# Crate matrices lists
create_matrices_list <- function(n, a, b) {
  matrix_list <- vector("list", n)
  for (i in 1:n) {
    matrix_list[[i]] <- matrix(0, ncol = a, nrow = b)
  }
  return(matrix_list)
}


