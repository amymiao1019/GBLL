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



# Calculate the MAPE for each horizon month h and each age group
mape_h <- function(n_pop, n_age, error_list, h_months, week_indicator){
  mape_list <- vector("list", n_pop)
  for(i in 1:n_pop){
    individual_error <- c()
    for(h in 1:h_months){
      mape <- colSums(error_list[[i]][1:week_indicator[h+1],])/week_indicator[h+1]
      individual_error <- rbind(individual_error, mape)
    }
  mape_list[[i]] <- individual_error
  }
  return(mape_list)
}



# Calculate the mean error for each country
mean_error_country <- function(error_list){
  country_split <- lapply(1:30, function(country_idx) {
    lapply(error_list, function(win) win[[country_idx]])
  })
  
  mean_errors_by_country <- lapply(country_split, function(mats) {
    mat_array <- simplify2array(mats)
    matrix(apply(mat_array, 1, mean), ncol = 1)
  })
  return(mean_errors_by_country)
}



# Calculate the MASE eror measure
mean_mase <- function(forecast_list, historical, week_indicator, horizon) {
  n_windows <- length(forecast_list)        # e.g. 10 expanding windows
  n_pop <- length(forecast_list[[1]])       # e.g. 30 countries
  n_age <- ncol(forecast_list[[1]][[1]])    # e.g. 4 age groups
  h <- nrow(forecast_list[[1]][[1]])        # e.g. 52 weeks
  
  n_months <- length(week_indicator) - 1
  mase_month_age <- matrix(0, nrow = n_months, ncol = n_age)
  count_valid <- matrix(0, nrow = n_months, ncol = n_age)
  
  # Loop over all windows and countries
  for (w in seq_len(n_windows)) {
    for (i in seq_len(n_pop)) {
      fc_mat <- forecast_list[[w]][[i]]
      actual_list <- lapply(historical, function(x) as.matrix(x[(training_period[w]+1):(training_period[w]+52),]))
      act_mat <- actual_list[[i]]
      hist_list <- lapply(historical, function(x) as.matrix(x[1:training_period[w],]))
      hist_mat <- hist_list[[i]]
      
      for (a in seq_len(n_age)) {
        hist_series <- hist_mat[, a]
        scale <- mean(abs(hist_series[(horizon + 1):length(hist_series)] -
                          hist_series[1:(length(hist_series) - horizon)]), na.rm = TRUE)
 
        # Loop over months
        for (h_idx in 1:n_months) {
          end_wk <- week_indicator[h_idx + 1]
          mae_h <- mean(abs(fc_mat[1:end_wk, a] -
                            act_mat[1:end_wk, a]), na.rm = TRUE)
          
          mase_month_age[h_idx, a] <- mase_month_age[h_idx, a] + (mae_h / scale)
          count_valid[h_idx, a] <- count_valid[h_idx, a] + 1
        }
      }
    }
  }
  
  # Average across countries × windows
  mase_month_age <- mase_month_age / count_valid
  rownames(mase_month_age) <- paste0("Month_", seq_len(n_months))
  colnames(mase_month_age) <- paste0("AgeGroup_", seq_len(n_age))
  
  return(mase_month_age)
}
