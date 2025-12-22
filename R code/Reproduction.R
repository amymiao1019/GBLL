1. Data preparation
date <- as.Date("2015-01-08") + 7*seq(0,52*5-1,1)
week_indicator <- c(0,4,9,13,17,22,26,30,35,39,43,48,52)
read_excel_name_v3 <- c("Australia", "Austria", "Belgium", "Bulgaria", "Canada", "Croatia", "Czech Republic", "Denmark", "England and Wales", "Finland", 
                        "France", "Germany", "Greece", "Hungary", "Italy", "Lithuania", "Netherlands", "New Zealand", "Norway", "Poland", "Portugal",
                        "Scotland", "Slovakia", "Slovenia", "South Korea", "Spain", "Sweden", "Switzerland", "Taiwan", "USA")

excel_files_v3 <- list.files("STMF dataset updated v3", pattern = "\\.xlsx$", full.names = TRUE)
data_list_v3 <- purrr::map(excel_files_v3, ~ readxl::read_excel(.x))
names(data_list_v3) <- basename(excel_files_v3)
name_HBY_v3 <- c("australia", "austria", "belgium", "bulgaria", "canada", "croatia", "czechrepublic", "denmark", "englandandwales", "finland", "france", 
                 "germany", "greece", "hungary", "italy", "lithuania", "netherlands", "newzealand", "norway", "poland", "portugal", "scotland", "slovakia", 
                 "slovenia", "southkorea", "spain", "sweden", "switzerland", "taiwan", "usa")
name_country_code_v3 <- c("AUS", "AUT", "BEL", "BGR", "CAN", "HRV", "CZE", "DNK", "ENW", "FIN", "FRA", "DEU", "GRC", "HUN", "ITA", "LTU", "NLD", "NZL",
                          "NOR", "POL", "PRT", "SCT", "SVK", "SVN", "KOR", "ESP", "SWE", "CHE", "TWN", "USA")
flip_indicator_whole_v3 <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
mr_list_v3 <- list()
for(i in 1:length(read_excel_name_v3)){
  matrix <- data_list_v3[[i]][,3:6]
  mr_list_v3[[i]] <- matrix
}



2. Do LL, HBY and GBLL for 10-fold expanding windows
## create lists to store the outputs for each expanding window
fitted_results_LL <- list()
error_LL <- list()
individual_error_LL <- list()
forecasted_mr_LL <- list()

fitted_results_HBY <- list()
error_HBY <- list()
individual_error_HBY <- list()
forecasted_mr_HBY <- list()

fitted_results_GBLL <- list()
error_GBLL <- list()
individual_error_GBLL <- list()
forecasted_mr_GBLL <- list()

## do 10 times of expanding window
training_period <- 208-week_indicator[1:10]
for(i in 1:10){
  # prepare the data needed
  train_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[1:training_period[i],]))
  train_mr_list_v3[[1]] <- 1/train_mr_list_v3[[1]]
  train_mr_list_v3[[18]] <- 1/train_mr_list_v3[[18]]
  train_mr_list_v3 <- lapply(train_mr_list_v3, log)
  test_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  train_regressor <- fourier_regressor(date, training_period[i], 52, 2)$train
  test_regressor <- fourier_regressor(date, training_period[i], 52, 2)$test
  
  # LL
  fitted_results_LL[[i]] <- gbll_original(length(read_excel_name_v3), 4, train_mr_list_v3, training_period[i],52,1)
  forecasted_LL <- gbll_forecasting(length(read_excel_name_v3), 4, 52, fitted_results_LL[[i]], train_regressor, test_regressor, flip_indicator_whole_v3, test_mr_list_v3, week_indicator, 12, fitted_results_LL[[i]]$iteration, 1)
  error_LL[[i]] <- forecasted_LL$error
  individual_error_LL[[i]] <- mape_h(length(read_excel_name_v3),4,forecasted_LL$individual_error,12,week_indicator)
  forecasted_mr_LL[[i]] <- forecasted_LL$forecasted_mortality
  
  # HBY
  fitted_results_HBY[[i]] <- hby_fit(name_HBY_v3, paste0("30 countries whole MR - EW", i-1, ".txt"), paste0("30 countries whole Exposure - EW", i-1, ".txt"), 6, length(read_excel_name_v3), 4, training_period[i])
  forecasted_HBY <- gbll_forecasting(length(read_excel_name_v3), 4, 52, fitted_results_HBY[[i]], train_regressor, test_regressor, flip_indicator_whole_v3, test_mr_list_v3, week_indicator, 12, 6, c(1,1,1,1,1,1))
  error_HBY[[i]] <- forecasted_HBY$error
  individual_error_HBY[[i]] <- mape_h(length(read_excel_name_v3),4,forecasted_HBY$individual_error,12,week_indicator)
  forecasted_mr_HBY[[i]] <- forecasted_HBY$forecasted_mortality
  
  #GBLL
  fitted_results_GBLL[[i]] <- gbll_original(length(read_excel_name_v3),4,train_mr_list_v3,training_period[i],52,50)
  forecasted_GBLL <- gbll_forecasting(length(read_excel_name_v3), 4, 52, fitted_results_GBLL[[i]], train_regressor, test_regressor, flip_indicator_whole_v3, test_mr_list_v3, week_indicator, 12, fitted_results_GBLL[[i]]$iteration, fitted_results_GBLL[[i]]$gamma)
  error_GBLL[[i]] <- forecasted_GBLL$error
  individual_error_GBLL[[i]] <- mape_h(length(read_excel_name_v3),4,forecasted_GBLL$individual_error,12,week_indicator)
  forecasted_mr_GBLL[[i]] <- forecasted_GBLL$forecasted_mortality
}

# Final mean error
mean_error_LL <- Reduce("+", unlist(error_LL, recursive = FALSE)) / (10*30)
mean_error_HBY <- Reduce("+", unlist(error_HBY, recursive = FALSE)) / (10*30)
mean_error_GBLL <- Reduce("+", unlist(error_GBLL, recursive = FALSE)) / (10*30)

# Calculate the mean MAPE for each age group
## LL
all_errors_LL <- unlist(individual_error_LL, recursive = FALSE)
combined_array_LL <- simplify2array(all_errors_LL)
mean_error_AG_LL <- apply(combined_array_LL, c(1, 2), mean)

## HBY
all_errors_HBY <- unlist(individual_error_HBY, recursive = FALSE)
combined_array_HBY <- simplify2array(all_errors_HBY)
mean_error_AG_HBY <- apply(combined_array_HBY, c(1, 2), mean)

## GBLL
all_errors_GBLL <- unlist(individual_error_GBLL, recursive = FALSE)
combined_array_GBLL <- simplify2array(all_errors_GBLL)
mean_error_AG_GBLL <- apply(combined_array_GBLL, c(1, 2), mean)




3. Prediction interval for the LL, HBY and GBLL model based on 100 simulations
error_LL_pi_lower <- list()
error_LL_pi_upper <- list()
forecasted_mr_LL_pi_lower <- list()
forecasted_mr_LL_pi_upper <- list()

error_HBY_pi_lower <- list()
error_HBY_pi_upper <- list()
forecasted_mr_HBY_pi_lower <- list()
forecasted_mr_HBY_pi_upper <- list()

error_GBLL_pi_lower <- list()
error_GBLL_pi_upper <- list()
forecasted_mr_GBLL_pi_lower <- list()
forecasted_mr_GBLL_pi_upper <- list()

# LL
for(i in 1:10){
  test_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  train_regressor <- fourier_regressor(date, training_period[i], 52, 2)$train
  test_regressor <- fourier_regressor(date, training_period[i], 52, 2)$test
  
  LL_pi_100 <- simulation_pi(30, 4, 52, fitted_results_LL[[i]], train_regressor, test_regressor,
                             flip_indicator_whole_v3, simulation_iter = 100, fitted_results_LL[[i]]$gamma, 
                             quantile = c(0.025, 0.975), 1019)
  forecasted_mr_LL_pi_lower[[i]] <- LL_pi_100$lower
  forecasted_mr_LL_pi_upper[[i]] <- LL_pi_100$upper
  
  error_LL_pi_lower[[i]] <- compute_forecast_error(forecasted_mr_LL_pi_lower[[i]], test_mr_list_v3,30,4,12,week_indicator)
  error_LL_pi_upper[[i]] <- compute_forecast_error(forecasted_mr_LL_pi_upper[[i]], test_mr_list_v3,30,4,12,week_indicator)
}

# HBY
for(i in 1:10){
  test_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  train_regressor <- fourier_regressor(date, training_period[i], 52, 2)$train
  test_regressor <- fourier_regressor(date, training_period[i], 52, 2)$test
  
  HBY_pi_100 <- simulation_pi(30, 4, 52, fitted_results_HBY[[i]], train_regressor, test_regressor,
                             flip_indicator_whole_v3, simulation_iter = 100, c(1,1,1,1,1,1), 
                             quantile = c(0.025, 0.975), 1019)
  forecasted_mr_HBY_pi_lower[[i]] <- HBY_pi_100$lower
  forecasted_mr_HBY_pi_upper[[i]] <- HBY_pi_100$upper
  
  error_HBY_pi_lower[[i]] <- compute_forecast_error(forecasted_mr_HBY_pi_lower[[i]], test_mr_list_v3,30,4,12,week_indicator)
  error_HBY_pi_upper[[i]] <- compute_forecast_error(forecasted_mr_HBY_pi_upper[[i]], test_mr_list_v3,30,4,12,week_indicator)
}

# GBLL
for(i in 1:10){
  test_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  train_regressor <- fourier_regressor(date, training_period[i], 52, 2)$train
  test_regressor <- fourier_regressor(date, training_period[i], 52, 2)$test
  
  GBLL_pi_100 <- simulation_pi(30, 4, 52, fitted_results_GBLL[[i]], train_regressor, test_regressor,
                             flip_indicator_whole_v3, simulation_iter = 100, fitted_results_GBLL[[i]]$gamma, 
                             quantile = c(0.025, 0.975), 1019)
  forecasted_mr_GBLL_pi_lower[[i]] <- GBLL_pi_100$lower
  forecasted_mr_GBLL_pi_upper[[i]] <- GBLL_pi_100$upper
  
  error_GBLL_pi_lower[[i]] <- compute_forecast_error(forecasted_mr_GBLL_pi_lower[[i]], test_mr_list_v3,30,4,12,week_indicator)
  error_GBLL_pi_upper[[i]] <- compute_forecast_error(forecasted_mr_GBLL_pi_upper[[i]], test_mr_list_v3,30,4,12,week_indicator)
}




4. To calculate the coverage based on the prediction interval
total_inside_LL <- 0
total_inside_HBY <- 0
total_inside_GBLL <- 0

# LL
for(i in 1:10){
  test_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  for(j in 1:30){
    lower <- forecasted_mr_LL_pi_lower[[i]][[j]]
    upper <- forecasted_mr_LL_pi_upper[[i]][[j]]
    actual <- test_mr_list_v3[[j]]
    inside <- (actual >= lower) & (actual <= upper)
    total_inside_LL <- total_inside_LL + sum(inside, na.rm = TRUE)
  }
}
total_inside_LL <- total_inside_LL/(52*4*30*10)

# HBY
for(i in 1:10){
  test_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  for(j in 1:30){
    lower <- forecasted_mr_HBY_pi_lower[[i]][[j]]
    upper <- forecasted_mr_HBY_pi_upper[[i]][[j]]
    actual <- test_mr_list_v3[[j]]
    inside <- (actual >= lower) & (actual <= upper)
    total_inside_HBY <- total_inside_HBY + sum(inside, na.rm = TRUE)
  }
}
total_inside_HBY <- total_inside_HBY/(52*4*30*10)

# GBLL
for(i in 1:10){
  test_mr_list_v3 <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  for(j in 1:30){
    lower <- forecasted_mr_GBLL_pi_lower[[i]][[j]]
    upper <- forecasted_mr_GBLL_pi_upper[[i]][[j]]
    actual <- test_mr_list_v3[[j]]
    inside <- (actual >= lower) & (actual <= upper)
    total_inside_GBLL <- total_inside_GBLL + sum(inside, na.rm = TRUE)
  }
}
total_inside_GBLL <- total_inside_GBLL/(52*4*30*10)

# summary
df_coverage <- data.frame(Model = c("LL", "HBY", "GBLL"),
                          Total_Inside = c(total_inside_LL, total_inside_HBY, total_inside_GBLL))




5. Clustering Method 1 results
# create lists to store the outputs for each expanding window
fitted_results_LL_m1_combined <- list()
error_LL_m1_combined <- list()
individual_error_LL_m1_combined <- list()

fitted_results_HBY_m1_combined <- list()
error_HBY_m1_combined <- list()
individual_error_HBY_m1_combined <- list()

fitted_results_GBLL_m1_combined <- list()
error_GBLL_m1_combined <- list()
individual_error_GBLL_m1_combined <- list()

# do 10 times of expanding window for each of the three clusters
for(j in 1:3){
  train_mr_list_v3_m1 <- list()
  test_mr_list_v3_m1 <- list()
  
  # create lists to store the outputs for each expanding window
  fitted_results_LL_m1 <- list()
  error_LL_m1 <- list()
  individual_error_LL_m1 <- list()

  fitted_results_HBY_m1 <- list()
  error_HBY_m1 <- list()
  individual_error_HBY_m1 <- list()

  fitted_results_GBLL_m1 <- list()
  error_GBLL_m1 <- list()
  individual_error_GBLL_m1 <- list()
  
  for(i in 1:10){
    index <- as.numeric(which(ordered_raw_lc_kapa_ts_NH_v3 == j))
    # prepare the data needed
    for(k in 1:length(index)){
      matrix <- as.matrix(data_list_v3[[index[k]]][1:training_period[i],3:6])
      train_mr_list_v3_m1[[k]] <- matrix
      test_mr_list_v3_m1[[k]] <- as.matrix(data_list_v3[[index[k]]][(training_period[i]+1):(training_period[i]+52),3:6])
    }
    
    flip_indicator <- rep(0, length(index))
    train_mr_list_v3_m1 <- lapply(train_mr_list_v3_m1, log)
    train_regressor <- fourier_regressor(date, training_period[i], 52, 2)$train
    test_regressor <- fourier_regressor(date, training_period[i], 52, 2)$test
    
    # LL
    fitted_results_LL_m1[[i]] <- gbll_original(length(index), 4, train_mr_list_v3_m1, training_period[i],52,1)
    forecasted_LL_m1 <- gbll_forecasting(length(index), 4, 52, fitted_results_LL_m1[[i]], train_regressor, test_regressor, flip_indicator, test_mr_list_v3_m1, week_indicator, 12, fitted_results_LL_m1[[i]]$iteration, 1)
    error_LL_m1[[i]] <- forecasted_LL_m1$error
    individual_error_LL_m1[[i]] <- mape_h(length(index),4,forecasted_LL_m1$individual_error,12,week_indicator)
    
    # HBY
    fitted_results_HBY_m1[[i]] <- hby_fit(name_HBY_v3[index], paste0("30 countries Clustering M5 C", j, " MR - EW", i - 1, ".txt"), paste0("30 countries Clustering M5 C", j, " Exposure - EW", i - 1, ".txt"), 6, length(index), 4, training_period[i])
    forecasted_HBY_m1 <- gbll_forecasting(length(index), 4, 52, fitted_results_HBY_m1[[i]], train_regressor, test_regressor, flip_indicator, test_mr_list_v3_m1, week_indicator, 12, 6, c(1,1,1,1,1,1))
    error_HBY_m1[[i]] <- forecasted_HBY_m1$error
    individual_error_HBY_m1[[i]] <- mape_h(length(index),4,forecasted_HBY_m1$individual_error,12,week_indicator)
    
    # GBLL
    fitted_results_GBLL_m1[[i]] <- gbll_original(length(index),4,train_mr_list_v3_m1,training_period[i],52,50)
    forecasted_GBLL_m1 <- gbll_forecasting(length(index), 4, 52, fitted_results_GBLL_m1[[i]], train_regressor, test_regressor, flip_indicator, test_mr_list_v3_m1, week_indicator, 12, fitted_results_GBLL_m1[[i]]$iteration, fitted_results_GBLL_m1[[i]]$gamma)
    error_GBLL_m1[[i]] <- forecasted_GBLL_m1$error
    individual_error_GBLL_m1[[i]] <- mape_h(length(index),4,forecasted_GBLL_m1$individual_error,12,week_indicator)
  }
  
  # combine the results across 3 clusters
  fitted_results_LL_m1_combined[[j]] <- fitted_results_LL_m1
  error_LL_m1_combined[[j]] <- error_LL_m1
  individual_error_LL_m1_combined[[j]] <- individual_error_LL_m1
  
  fitted_results_HBY_m1_combined[[j]] <- fitted_results_HBY_m1
  error_HBY_m1_combined[[j]] <- error_HBY_m1
  individual_error_HBY_m1_combined[[j]] <- individual_error_HBY_m1
  
  fitted_results_GBLL_m1_combined[[j]] <- fitted_results_GBLL_m1
  error_GBLL_m1_combined[[j]] <- error_GBLL_m1
  individual_error_GBLL_m1_combined[[j]] <- individual_error_GBLL_m1
}


# final mean error
mean_error_LL_m1 <- apply(simplify2array(flatten(flatten(error_LL_m1_combined))), c(1, 2), mean)
mean_error_HBY_m1 <- apply(simplify2array(flatten(flatten(error_HBY_m1_combined))), c(1, 2), mean)
mean_error_GBLL_m1 <- apply(simplify2array(flatten(flatten(error_GBLL_m1_combined))), c(1, 2), mean)






6. Clustering Method 2 results
# create lists to store the outputs for each expanding window
fitted_results_LL_m2_combined <- list()
error_LL_m2_combined <- list()
individual_error_LL_m2_combined <- list()

fitted_results_HBY_m2_combined <- list()
error_HBY_m2_combined <- list()
individual_error_HBY_m2_combined <- list()

fitted_results_GBLL_m2_combined <- list()
error_GBLL_m2_combined <- list()
individual_error_GBLL_m2_combined <- list()

# do 10 times of expanding window for each of the three clusters
for(j in 1:3){
  train_mr_list_v3_m2 <- list()
  test_mr_list_v3_m2 <- list()
  
  # create lists to store the outputs for each expanding window
  fitted_results_LL_m2 <- list()
  error_LL_m2 <- list()
  individual_error_LL_m2 <- list()

  fitted_results_HBY_m2 <- list()
  error_HBY_m2 <- list()
  individual_error_HBY_m2 <- list()

  fitted_results_GBLL_m2 <- list()
  error_GBLL_m2 <- list()
  individual_error_GBLL_m2 <- list()
  
  for(i in 1:10){
    index <- as.numeric(which(lc_kapa_stl_trend_v3$cluster == j))
    # prepare the data needed
    for(k in 1:length(index)){
      matrix <- as.matrix(data_list_v3[[index[k]]][1:training_period[i],3:6])
      train_mr_list_v3_m2[[k]] <- matrix
      test_mr_list_v3_m2[[k]] <- as.matrix(data_list_v3[[index[k]]][(training_period[i]+1):(training_period[i]+52),3:6])
    }
    
    flip_indicator <- rep(0, length(index))
    flip_indicator[index == 1 | index == 18] <- 1
    which_flip <- which(flip_indicator==1)
    if(length(which_flip)==0){
    } else{
        train_mr_list_v3_m2[[which_flip]] <- 1/train_mr_list_v3_m2[[which_flip]]
      }
    train_mr_list_v3_m2 <- lapply(train_mr_list_v3_m2, log)
    train_regressor <- fourier_regressor(date, training_period[i], 52, 2)$train
    test_regressor <- fourier_regressor(date, training_period[i], 52, 2)$test
    
    # LL
    fitted_results_LL_m2[[i]] <- gbll_original(length(index), 4, train_mr_list_v3_m2, training_period[i],52,1)
    forecasted_LL_m2 <- gbll_forecasting(length(index), 4, 52, fitted_results_LL_m2[[i]], train_regressor, test_regressor, flip_indicator, test_mr_list_v3_m2, week_indicator, 12, fitted_results_LL_m2[[i]]$iteration, 1)
    error_LL_m2[[i]] <- forecasted_LL_m2$error
    individual_error_LL_m2[[i]] <- mape_h(length(index),4,forecasted_LL_m2$individual_error,12,week_indicator)
    
    # HBY
    fitted_results_HBY_m2[[i]] <- hby_fit(name_HBY_v3[index], paste0("30 countries Clustering M2 C", j, " MR - EW", i - 1, ".txt"), paste0("30 countries Clustering M2 C", j, " Exposure - EW", i - 1, ".txt"), 6, length(index), 4, training_period[i])
    forecasted_HBY_m2 <- gbll_forecasting(length(index), 4, 52, fitted_results_HBY_m2[[i]], train_regressor, test_regressor, flip_indicator, test_mr_list_v3_m2, week_indicator, 12, 6, c(1,1,1,1,1,1))
    error_HBY_m2[[i]] <- forecasted_HBY_m2$error
    individual_error_HBY_m2[[i]] <- mape_h(length(index),4,forecasted_HBY_m2$individual_error,12,week_indicator)
    
    # GBLL
    fitted_results_GBLL_m2[[i]] <- gbll_original(length(index),4,train_mr_list_v3_m2,training_period[i],52,50)
    forecasted_GBLL_m2 <- gbll_forecasting(length(index), 4, 52, fitted_results_GBLL_m2[[i]], train_regressor, test_regressor, flip_indicator, test_mr_list_v3_m2, week_indicator, 12, fitted_results_GBLL_m2[[i]]$iteration, fitted_results_GBLL_m2[[i]]$gamma)
    error_GBLL_m2[[i]] <- forecasted_GBLL_m2$error
    individual_error_GBLL_m2[[i]] <- mape_h(length(index),4,forecasted_GBLL_m2$individual_error,12,week_indicator)
  }
  
  # combine the results across 3 clusters
  fitted_results_LL_m2_combined[[j]] <- fitted_results_LL_m2
  error_LL_m2_combined[[j]] <- error_LL_m2
  individual_error_LL_m2_combined[[j]] <- individual_error_LL_m2
  
  fitted_results_HBY_m2_combined[[j]] <- fitted_results_HBY_m2
  error_HBY_m2_combined[[j]] <- error_HBY_m2
  individual_error_HBY_m2_combined[[j]] <- individual_error_HBY_m2
  
  fitted_results_GBLL_m2_combined[[j]] <- fitted_results_GBLL_m2
  error_GBLL_m2_combined[[j]] <- error_GBLL_m2
  individual_error_GBLL_m2_combined[[j]] <- individual_error_GBLL_m2
}


# final mean error
mean_error_LL_m2 <- apply(simplify2array(flatten(flatten(error_LL_m2_combined))), c(1, 2), mean)
mean_error_HBY_m2 <- apply(simplify2array(flatten(flatten(error_HBY_m2_combined))), c(1, 2), mean)
mean_error_GBLL_m2 <- apply(simplify2array(flatten(flatten(error_GBLL_m2_combined))), c(1, 2), mean)
