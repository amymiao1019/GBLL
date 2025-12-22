# Fit the single-population Vector Autoregressive models to country-specific mortality rates
library(vars)
all_mape <- vector("list", 10) 
all_forecasts <- vector("list", 10)

for(i in 1:10){
    train_list <- lapply(mr_list_v3, function(x) as.matrix(x[1:training_period[i], ]))
    window_forecasts <- vector("list", length(train_list))
    
    for (c in seq_len(length(train_list))) {
        y <- train_list[[c]]
        y_log <- log(y)
        y_ts <- ts(y_log, frequency = 52)

        # stationarity check
        d_vec <- sapply(1:4, function(k) ndiffs(na.omit(y_ts[, k])))
        D_vec <- sapply(1:4, function(k) nsdiffs(na.omit(y_ts[, k])))
        d <- min(max(d_vec),1)
        D <- min(max(D_vec),1)
        
        # apply differencing
        y_diff <- y_ts
        if (D > 0) y_diff <- diff(y_diff, lag = 52, differences = D)
        if (d > 0) y_diff <- diff(y_diff, differences = d)
        y_diff <- na.omit(y_diff)
        
        # lag selection
        lag_selection <- VARselect(y_diff, lag.max = 10, type = "const")
        p_opt <- lag_selection$selection["HQ(n)"]
        var_model <- VAR(y_diff, p = p_opt, type = "const")
        
        
        # forecast 52 weeks ahead
        fc <- predict(var_model, n.ahead = 52)
        fc_mat <- sapply(fc$fcst, function(x) x[, "fcst"])
        
        # revert back to original mortality
        last_log <- log(y[nrow(y), ])
        fc_level <- matrix(NA, nrow = 52, ncol = dim(train_list[[c]])[2])
        for (age in 1:dim(train_list[[c]])[2]) {
          if (d == 0 & D == 0) {
            # 1) no differencing
            fc_level[, age] <- exp(fc_mat[, age])
          } else if (d == 1 & D == 0) {
            # 2) regular difference only
            fc_level[, age] <- exp(last_log[age] + cumsum(fc_mat[, age]))
          } else if (d == 0 & D == 1) {
            # 3) seasonal difference only
            hist_window <- log(y[(nrow(y)-51):nrow(y), age])  # last 52 weeks
            log_preds <- numeric(52)
            for (t in 1:52) {
                log_preds[t] <- hist_window[t] + fc_mat[t, age]
              }
            fc_level[, age] <- exp(log_preds)
          } else if (d == 1 & D == 1) {
            # 4) both regular and seasonal differencing
            hist_window <- log(y[(nrow(y)-52):nrow(y), age]) # last 53 weeks
            y_season_diff <- numeric(52)
            y_forecast <- numeric(52)
            for (t in 1:52) {
              if (t == 1) {
                y_season_diff[t] <- hist_window[53] + fc_mat[t, age]  # last observed seasonal-diff
                y_forecast[t] <- y_season_diff[t] + hist_window[2] - hist_window[1]
              } else {
                y_season_diff[t] <- y_season_diff[t-1] + fc_mat[t, age]
                y_forecast[t] <- y_season_diff[t] + hist_window[t+1] - hist_window[t]
              }
            }
            fc_level[, age] <- exp(y_forecast)
          }
        }
        
        window_forecasts[[c]] <- fc_level
    }
    all_forecasts[[i]] <- window_forecasts
    
    # calculating the mape
    forecast_error_individual_var <- create_matrices_list(length(train_list), dim(train_list[[c]])[2], 52)
    rowsum_forecast_error_var <- create_matrices_list(length(train_list), 1, 52)
    forecast_error_var <- create_matrices_list(length(train_list), 1, 12)
    
    test_mortality <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
    difference_var <- Map(`-`, all_forecasts[[i]], test_mortality)
    abs_difference_var <- lapply(difference_var, function(diff) abs(diff))
    forecast_error_individual_var <- Map(`/`, abs_difference_var, test_mortality)
    
    rowsum_forecast_error_var <- lapply(forecast_error_individual_var, function(error) rowSums(error))
    
    for(m in 1:length(train_list)){
      for(n in 1:12){
        forecast_error_var[[m]][n] <- sum(rowsum_forecast_error_var[[m]][1:week_indicator[n+1]])/(dim(train_list[[c]])[2]*week_indicator[n+1])
      }
    }
    
    error_array <- array(unlist(forecast_error_var), dim = c(12, 1, length(forecast_error_var)))
    mean_error <- apply(error_array, 1, mean)
    mean_error <- matrix(mean_error, nrow = 12, ncol = 1)
    all_mape[[i]] <- mean_error
}

# find the mean
mean_mape_var <- Reduce("+", all_mape)/length(all_mape)
