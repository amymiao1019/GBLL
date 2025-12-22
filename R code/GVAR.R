# fit and forecast using the multi-population Global VAR model for raw mortality rates
library(GVARX)

# Pre-allocate
all_mape_gvar <- vector("list", 10)
all_forecasts_gvar <- vector("list", 10)

for(i in 1:10){
  train_list <- lapply(mr_list_v3, function(x) as.matrix(x[1:training_period[i], ]))
  train_list[[1]] <- 1 / train_list[[1]]
  train_list[[18]] <- 1 / train_list[[18]]
  train_list_log <- lapply(train_list, log)
  
  d_common <- 1
  D_common <- 1
  train_list_diff <- lapply(train_list_log, function(y_ts){
    y_diff <- diff(y_ts, differences = d_common)
    y_diff <- diff(y_diff, lag = 52, differences = D_common)
    return(na.omit(y_diff))
  })
  
  # create panel data
  panel_df <- data.frame()

  N <- length(train_list_diff)
  T_diff <- nrow(train_list_diff[[1]])
  K <- ncol(train_list_diff[[1]])
  
  for(c in 1:N){
    df_c <- as.data.frame(train_list_diff[[c]])
    df_c$Time <- date[54:training_period[i]]
    df_c$ID <- c  # country index
    panel_df <- rbind(panel_df, df_c)
  }
  panel_df <- panel_df[, c("ID", "Time", colnames(df_c)[1:K])]

  # Construct the weight
  W <- matrix(1/(N-1), nrow = N, ncol = N)
  diag(W) <- 0
  
  # fit the Global VAR model
  gvar_model <- GVARest(data = panel_df, lag.max = 10, type = "const", ic="HQ", weight.matrix = W)
  gvar_GF <- GVAR_GF(data = panel_df, p=gvar_model$p, type = "const", ic="HQ", weight.matrix = W)
  
  # forecast
  K_total <- nrow(last_obs)
  a0 <- numeric(K_total)
  p_max <- if(!is.null(gvar_GF$G2)) 2 else 1
  for(c in 1:N){
    coefs_c <- coef(gvar_model$gvar[[c]])  # list of matrices, one per domestic variable
    
    for(age_idx in 1:K){
      var_name <- names(coefs_c)[age_idx]       # e.g., "X1.75.84"
      const_value <- coefs_c[[var_name]]["const", "Estimate"]
      
      # Position in stacked vector
      row_index <- (c-1)*K + age_idx
      a0[row_index] <- const_value
    }
  }
  
  last_obs <- matrix(NA, nrow=K_total, ncol=p_max)
  for(c in 1:N){
    for(age in 1:K){
      row_index <- (c-1)*K + age
      last_obs[row_index, 1] <- train_list_diff[[c]][T_diff, age]       # X_t
      if(p_max==2){
        last_obs[row_index, 2] <- train_list_diff[[c]][T_diff-1, age]  # X_{t-1}
      }
    }
  }
  
  G0 <- gvar_GF$G0
  G1 <- gvar_GF$G1
  F1 <- gvar_GF$F1
  G2 <- if(!is.null(gvar_GF$G2)) gvar_GF$G2 else matrix(0, K_total, K_total)
  F2 <- if(!is.null(gvar_GF$F2)) gvar_GF$F2 else matrix(0, K_total, K_total)
  horizon <- 52
  seasonal_lag <- 52
  countries_inverted <- c(1,18)  # AUS/NZL
  
  G0_inv <- solve(G0)
  
  
  ## actual forecasting
  X_forecast <- matrix(NA, nrow=K_total, ncol=horizon)
  X_t   <- last_obs[,1]
  X_t_1 <- if(p_max==2) last_obs[,2] else rep(0, K_total)
  
  for(h in 1:horizon){
    # Foreign lags (weighted sum of other countries)
    X_foreign_lag1 <- kronecker(W, diag(K)) %*% X_t
    X_foreign_lag2 <- kronecker(W, diag(K)) %*% X_t_1
    
    # Forecast next step
    X_next <- G0_inv %*% a0 +
              G0_inv %*% G1 %*% X_t +
              G0_inv %*% G2 %*% X_t_1
    
    # Store forecast
    X_forecast[,h] <- X_next
    
    # Update lags
    if(p_max==2){
      X_t_1 <- X_t
    }
    X_t <- X_next
  }
  
  # convert back to original mortality rate
  forecast_list <- vector("list", N)
  for (c in 1:N) {
    y <- train_list[[c]]          
    ylog <- train_list_log[[c]]   
    T0 <- nrow(ylog)
    
    fc_mat <- matrix(NA, nrow = horizon, ncol = K)
    
    for (age in 1:K) {
      row_index <- (c-1)*K + age
      fc_mat[, age] <- X_forecast[row_index, ]   # model output for this variable
    }
    
    fc_level <- matrix(NA, nrow = horizon, ncol = K)
    for (age in 1:K) {
      # 53-week history window in log-scale
      hist_window <- ylog[(T0-horizon):T0, age]   # length 53 (index 1...53)
      
      y_season_diff <- numeric(horizon)  
      y_forecast <- numeric(horizon)
      
      for (t in 1:horizon) {
        if (t == 1) {
          y_season_diff[t] <- hist_window[53] + fc_mat[t, age]
          y_forecast[t] <- y_season_diff[t] + hist_window[2] - hist_window[1]
        } else {
          y_season_diff[t] <- y_season_diff[t-1] + fc_mat[t, age]
          y_forecast[t] <- y_season_diff[t] + hist_window[t+1] - hist_window[t]
        }
      }
      
      if (c %in% countries_inverted) {
        fc_level[, age] <- 1 / exp(y_forecast)
      } else {
        fc_level[, age] <- exp(y_forecast)
      }
    }
    forecast_list[[c]] <- fc_level
  }
  
  all_forecasts_gvar[[i]] <- forecast_list
    
  # calculating the mape
  forecast_error_individual_gvar <- create_matrices_list(length(train_list), dim(train_list[[c]])[2], 52)
  rowsum_forecast_error_gvar <- create_matrices_list(length(train_list), 1, 52)
  forecast_error_gvar <- create_matrices_list(length(train_list), 1, 12)
    
  test_mortality <- lapply(mr_list_v3, function(x) as.matrix(x[(training_period[i]+1):(training_period[i]+52),]))
  difference_gvar <- Map(`-`, all_forecasts_gvar[[i]], test_mortality)
  abs_difference_gvar <- lapply(difference_gvar, function(diff) abs(diff))
  forecast_error_individual_gvar <- Map(`/`, abs_difference_gvar, test_mortality)
    
  rowsum_forecast_error_gvar <- lapply(forecast_error_individual_gvar, function(error) rowSums(error))
    
  for(m in 1:length(train_list)){
    for(n in 1:12){
      forecast_error_gvar[[m]][n] <- sum(rowsum_forecast_error_gvar[[m]][1:week_indicator[n+1]])/(dim(train_list[[c]])[2]*week_indicator[n+1])
    }
  }
    
  error_array_gvar <- array(unlist(forecast_error_gvar), dim = c(12, 1, length(forecast_error_gvar)))
  mean_error_gvar <- apply(error_array_gvar, 1, mean)
  mean_error_gvar <- matrix(mean_error_gvar, nrow = 12, ncol = 1)
  all_mape_gvar[[i]] <- mean_error_gvar
}

mean_mape_gvar <- Reduce("+", all_mape_gvar)/length(all_mape_gvar)
