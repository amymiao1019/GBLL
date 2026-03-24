1. SARIMA
gbll_forecasting_SARIMA <- function(n_pop, n_age, horizon, fitted_output, flip_indicator, test_mortality, week_indicator, h_months, iteration, gamma){
  
  # rename and store the values of ax, bx, kt
  ax <- fitted_output$ax
  bx <- fitted_output$bx
  kt <- fitted_output$kapa
  
  kt_forecast <- create_matrices_list((n_pop+1), iteration, horizon)
  
  # forecast the kt
  for(d in 1:(n_pop+1)){
    for(f in 1:iteration){
      fit.kt <- auto.arima(ts(kt[[d]][,f], frequency = horizon), seasonal = TRUE, max.p = 2, max.q = 2, max.P = 1, max.Q = 1, max.d=1, max.D=1, stepwise = TRUE, approximation = TRUE)
      kt_forecast[[d]][,f] <- forecast(fit.kt, h = horizon)$mean
    }
  }
  
  # find the log forecasted mortality
  forecasted_mortality <- create_matrices_list(n_pop, n_age, horizon)
  for(z in 1:iteration){
    for(g in 1:n_pop){
      for(h in 1:n_age){
        forecasted_mortality[[g]][,h] <- forecasted_mortality[[g]][,h] + gamma[z]*(ax[[1]][z,h] + bx[[1]][z,h]*kt_forecast[[1]][,z] + ax[[g+1]][z,h] + bx[[g+1]][z,h]*kt_forecast[[g+1]][,z])
      }
    }
  }
 
  
  # find the actual forecasted mortality
  for(p in 1:n_pop){
    if(flip_indicator[p]==0){
      forecasted_mortality[[p]] <- exp(forecasted_mortality[[p]])
    } else{
      forecasted_mortality[[p]] <- 1/exp(forecasted_mortality[[p]])
    }
  }
  
  
  # find the forecasting error
  forecast_error_individual <- create_matrices_list(n_pop, n_age, horizon)
  rowsum_forecast_error <- create_matrices_list(n_pop, 1, horizon)
  rowsum_mae <- create_matrices_list(n_pop, 1, horizon)
  forecast_error <- create_matrices_list(n_pop, 1, h_months)
  forecast_mae <- create_matrices_list(n_pop, 1, h_months)
  
  difference <- Map(`-`, forecasted_mortality, test_mortality)
  abs_difference <- lapply(difference, function(diff) abs(diff))
  forecast_error_individual <- Map(`/`, abs_difference, test_mortality)
  
  rowsum_forecast_error <- lapply(forecast_error_individual, function(error) rowSums(error))
  rowsum_mae <- lapply(abs_difference, function(error) rowSums(error))
  
  for(m in 1:n_pop){
    for(n in 1:h_months){
      forecast_error[[m]][n] <- sum(rowsum_forecast_error[[m]][1:week_indicator[n+1]])/(n_age*week_indicator[n+1])
      forecast_mae[[m]][n] <- sum(rowsum_mae[[m]][1:week_indicator[n+1]])/(n_age*week_indicator[n+1])
    }
  }
  
  
  # return the value
  output <- list(forecasted_mortality=forecasted_mortality, error=forecast_error, individual_error = forecast_error_individual, mae=forecast_mae)
  return(output)
}





2. Absolute loss function
# loss function
obj_abs <- function(n_pop, n_age, mortality_rate, fitted_mortality, gamma){
  obj_sum <- 0
  
  if(n_pop == 1){
    for(j in 1:n_age){
      fitted <- fitted_mortality[, j]
      error <- mortality_rate[, j] - gamma * fitted
      abs_error <- sum(abs(error))    # use L1 loss
      obj_sum <- obj_sum + abs_error
    }
    
  } else {
    for(i in 1:n_pop){
      for(j in 1:n_age){
        fitted <- fitted_mortality[[i]][, j]
        error <- mortality_rate[[i]][, j] - gamma * fitted
        abs_error <- sum(abs(error))  # use L1 loss
        obj_sum <- obj_sum + abs_error
      }
    }
  }
  
  return(obj_sum)
}





gbll_original_abs_loss_function <- function(n_pop, n_age, mortality_rate_list, train_time, horizon, max_iteration){
  # create the initial values
  log_residual <- mortality_rate_list
  iteration <- max_iteration
  gamma <- c()
  
  # create the ax, bx, kt, p_value, p_value indicator, drop indicator matrix list to store values
  ax <- create_matrices_list((n_pop+1), n_age, max_iteration)
  bx <- create_matrices_list((n_pop+1), n_age, max_iteration)
  kt <- create_matrices_list((n_pop+1), max_iteration, train_time)
  p_value <- create_matrices_list(n_pop, n_age, max_iteration)
  stop.indicator <- create_matrices_list(n_pop, n_age, max_iteration)
  fitted_value <- create_matrices_list(n_pop, n_age, train_time)
  
  
  for(i in 1:max_iteration){
    sum.p_value <- 0
    # set the product and ratio matrices
    exp_residual <- lapply(log_residual, function(mat) exp(mat))
    product <- (exp_residual |> reduce(`*`))^(1/n_pop)
    ratio_list <- mapply(function(mat_list) {mat_list/product}, exp_residual, SIMPLIFY = FALSE)
    ratio_list_log <- lapply(ratio_list, function(matrix) log(matrix))
    
    # now fit the LC on product and each ratio
    product_fitted <- lc(log(product))
    ratio_fitted <- lapply(ratio_list_log, lc)
    
    # update the new residual after LC and find Gamma
    new_log_residual <- list()
    fitted <- list()
    for(j in 1:n_pop){
      matrix <- product_fitted$fitted + ratio_fitted[[j]]$fitted
      fitted[[j]] <- matrix
    }
    
    gamma_value <- optim(par=1, obj_abs, method = "BFGS", n_pop = n_pop, n_age = n_age, mortality_rate = log_residual, fitted_mortality = fitted)$par
    gamma <- c(gamma, gamma_value)
    
    
    for(j in 1:n_pop){
      matrix_residual <- log_residual[[j]] - gamma_value* fitted[[j]]
      new_log_residual[[j]] <- matrix_residual
    }
    log_residual <- new_log_residual
    
    # store the values of ax, bx and kt
    ax[[1]][i,] <- product_fitted$ax
    bx[[1]][i,] <- product_fitted$bx
    kt[[1]][,i] <- product_fitted$kapa
    
    for(k in 1:n_pop){
      ax[[k+1]][i,] <- ratio_fitted[[k]]$ax
      bx[[k+1]][i,] <- ratio_fitted[[k]]$bx 
      kt[[k+1]][,i] <- ratio_fitted[[k]]$kapa
    }
    
    # calculate the p_value
    for(q in 1:n_pop){
      for(w in 1:n_age){
        p_value[[q]][i,w] <- Box.test(log_residual[[q]][,w], lag = horizon, type = "Ljung-Box")$p.value
        stop.indicator[[q]][i,w] <- ifelse(p_value[[q]][i,w] < 0.05, 0, 1)
      }
      sum.p_value <- sum.p_value + rowSums(stop.indicator[[q]])[i]
    }
    
    # update next iteration
    if(sum.p_value == n_pop * n_age){
      iteration <- i
      break
    }
  }
  
  # fitted values
  for(z in 1:iteration){
    for(g in 1:n_pop){
      for(h in 1:n_age){
        fitted_value[[g]][,h] <- fitted_value[[g]][,h] + gamma[z]*(ax[[1]][z,h] + bx[[1]][z,h]*kt[[1]][,z] + ax[[g+1]][z,h] + bx[[g+1]][z,h]*kt[[g+1]][,z])
      }
    }
  }
  
  for(j in 1:n_pop){
    fitted_value[[j]]=exp(fitted_value[[j]])
  }
  
  
  # return the output
  output <- list(ax=ax, kapa=kt, bx=bx, iteration=iteration, gamma=gamma, fitted_value=fitted_value)
  return(output)
}






3. Stopping criteria based on no improvement
gbll_original_stopping <- function(n_pop, n_age, mortality_rate_list, train_time,
                          horizon, max_iteration, epsilon = 1e-6) {
  
  # initial residuals
  log_residual <- mortality_rate_list
  gamma <- c()
  iteration <- max_iteration
  
  # storage
  ax <- create_matrices_list((n_pop+1), n_age, max_iteration)
  bx <- create_matrices_list((n_pop+1), n_age, max_iteration)
  kt <- create_matrices_list((n_pop+1), max_iteration, train_time)
  fitted_value <- create_matrices_list(n_pop, n_age, train_time)
  loss_vec <- numeric(max_iteration)
  
  for(i in 1:max_iteration){
    exp_residual <- lapply(log_residual, exp)
    product <- (exp_residual |> reduce(`*`))^(1/n_pop)
    ratio_list <- mapply(function(mat) mat / product, exp_residual, SIMPLIFY=FALSE)
    ratio_list_log <- lapply(ratio_list, log)
    
    product_fitted <- lc(log(product))
    ratio_fitted <- lapply(ratio_list_log, lc)
    
    fitted <- list()
    for(j in 1:n_pop){
      fitted[[j]] <- product_fitted$fitted + ratio_fitted[[j]]$fitted
    }
    
    gamma_value <- optim(
      par = 1,
      fn = obj,
      method = "BFGS",
      n_pop = n_pop,
      n_age = n_age,
      mortality_rate = log_residual,
      fitted_mortality = fitted
    )$par
    
    gamma <- c(gamma, gamma_value)
    
    log_residual <- lapply(1:n_pop, function(j) {
      log_residual[[j]] - gamma_value * fitted[[j]]
    })
    loss_vec[i] <- sum(unlist(log_residual)^2)

    
    ax[[1]][i,] <- product_fitted$ax
    bx[[1]][i,] <- product_fitted$bx
    kt[[1]][,i] <- product_fitted$kapa
    
    for(k in 1:n_pop){
      ax[[k+1]][i,] <- ratio_fitted[[k]]$ax
      bx[[k+1]][i,] <- ratio_fitted[[k]]$bx
      kt[[k+1]][,i] <- ratio_fitted[[k]]$kapa
    }
    
    # new stopping criteria
    if (i > 1) {
      if (abs(loss_vec[i] - loss_vec[i-1]) < epsilon) {
        iteration <- i
        break
      }
    }
  } # end loop
  
  
  for(z in 1:iteration){
    for(g in 1:n_pop){
      for(h in 1:n_age){
        fitted_value[[g]][,h] <- fitted_value[[g]][,h] +
          gamma[z] * (ax[[1]][z,h] + bx[[1]][z,h]*kt[[1]][,z] +
                        ax[[g+1]][z,h] + bx[[g+1]][z,h]*kt[[g+1]][,z])
      }
    }
  }
  
  fitted_value <- lapply(fitted_value, exp)
  return(list(
    ax = ax,
    kapa = kt,
    bx = bx,
    iteration = iteration,
    gamma = gamma,
    fitted_value = fitted_value,
    loss = loss_vec[1:iteration]
  ))
}





4. Stopping criteria based on AIC and/or BIC
gbll_likelihood <- function(n_pop, n_age, mortality_rate_list, train_time, horizon,
                                 max_iteration, criterion = c("AIC","BIC")) {
  
  criterion <- match.arg(criterion)
  
  log_residual <- mortality_rate_list
  iteration <- max_iteration
  gamma <- c()
  
  # IC tracking
  AICv <- rep(NA, max_iteration)
  BICv <- rep(NA, max_iteration)
  loglik_v <- rep(NA, max_iteration)
  RSSv <- rep(NA, max_iteration)
  
  # store ax, bx, kt
  ax <- create_matrices_list((n_pop+1), n_age, max_iteration)
  bx <- create_matrices_list((n_pop+1), n_age, max_iteration)
  kt <- create_matrices_list((n_pop+1), max_iteration, train_time)
  fitted_value <- create_matrices_list(n_pop, n_age, train_time)
  
  
  for(i in 1:max_iteration){
    exp_residual <- lapply(log_residual, function(mat) exp(mat))
    product <- (exp_residual |> reduce(`*`))^(1/n_pop)
    ratio_list <- mapply(function(mat_list) {mat_list/product}, exp_residual, SIMPLIFY = FALSE)
    ratio_list_log <- lapply(ratio_list, function(matrix) log(matrix))
    
    # now fit the LC on product and each ratio
    product_fitted <- lc(log(product))
    ratio_fitted <- lapply(ratio_list_log, lc)
    
    # update the new residual after LC and find Gamma
    new_log_residual <- list()
    fitted <- list()
    for(j in 1:n_pop){
      matrix <- product_fitted$fitted + ratio_fitted[[j]]$fitted
      fitted[[j]] <- matrix
    }
    
    gamma_value <- optim(par=1, obj, method = "BFGS", n_pop = n_pop, n_age = n_age, 
                         mortality_rate = log_residual, fitted_mortality = fitted)$par
    gamma <- c(gamma, gamma_value)
    
    
    for(j in 1:n_pop){
      matrix_residual <- log_residual[[j]] - gamma_value* fitted[[j]]
      new_log_residual[[j]] <- matrix_residual
    }
    log_residual <- new_log_residual
    
    # store the values of ax, bx and kt
    ax[[1]][i,] <- product_fitted$ax
    bx[[1]][i,] <- product_fitted$bx
    kt[[1]][,i] <- product_fitted$kapa
    
    for(k in 1:n_pop){
      ax[[k+1]][i,] <- ratio_fitted[[k]]$ax
      bx[[k+1]][i,] <- ratio_fitted[[k]]$bx 
      kt[[k+1]][,i] <- ratio_fitted[[k]]$kapa
    }
    
    # Compute AIC/BIC for iteration i
    resid_concat <- unlist(lapply(log_residual, function(m) as.numeric(m)))
    resid_concat <- resid_concat[!is.na(resid_concat)]
    RSS_i <- sum(resid_concat^2)
    Nobs <- length(resid_concat)
    n_param <- i*((n_pop + 1)*(2*n_age + train_time) + 1)
    
    sigma2_hat <- RSS_i / Nobs
    loglik_i <- - (Nobs/2) * (log(2*pi) + 1 + log(sigma2_hat))
    
    RSSv[i] <- RSS_i
    loglik_v[i] <- loglik_i
    AICv[i] <- -2 * loglik_i + 2 * n_param
    BICv[i] <- -2 * loglik_i + n_param * log(Nobs)
  }
  
  best_i_AIC <- which.min(AICv)
  best_i_BIC <- which.min(BICv)
  iteration <- ifelse(criterion == "AIC", best_i_AIC, best_i_BIC)  

  
  # fitted values
  for(z in 1:iteration){
    for(g in 1:n_pop){
      for(h in 1:n_age){
        fitted_value[[g]][,h] <- fitted_value[[g]][,h] +
          gamma[z] * ( ax[[1]][z,h] + bx[[1]][z,h]*kt[[1]][,z] +
                       ax[[g+1]][z,h] + bx[[g+1]][z,h]*kt[[g+1]][,z] )
      }
    }
  }
  
  for(j in 1:n_pop){
    fitted_value[[j]] <- exp(fitted_value[[j]])
  }
  
  return(list(
    ax=ax, kapa=kt, bx=bx,
    gamma=gamma,
    iteration=iteration,
    AIC_value = AICv,
    BIC_value = BICv,
    fitted_value=fitted_value
  ))
}

                         




5. Seasonal adjustment for Southern Hemisphere countries
gbll_original_stopping <- function(n_pop, n_age, mortality_rate_list, train_time,
                          horizon, max_iteration, epsilon = 1e-6) {
  
  # initial residuals
  log_residual <- mortality_rate_list
  gamma <- c()
  iteration <- max_iteration
  
  # storage
  ax <- create_matrices_list((n_pop+1), n_age, max_iteration)
  bx <- create_matrices_list((n_pop+1), n_age, max_iteration)
  kt <- create_matrices_list((n_pop+1), max_iteration, train_time)
  fitted_value <- create_matrices_list(n_pop, n_age, train_time)
  loss_vec <- numeric(max_iteration)
  
  for(i in 1:max_iteration){
    exp_residual <- lapply(log_residual, exp)
    product <- (exp_residual |> reduce(`*`))^(1/n_pop)
    ratio_list <- mapply(function(mat) mat / product, exp_residual, SIMPLIFY=FALSE)
    ratio_list_log <- lapply(ratio_list, log)
    
    product_fitted <- lc(log(product))
    ratio_fitted <- lapply(ratio_list_log, lc)
    
    fitted <- list()
    for(j in 1:n_pop){
      fitted[[j]] <- product_fitted$fitted + ratio_fitted[[j]]$fitted
    }
    
    gamma_value <- optim(
      par = 1,
      fn = obj,
      method = "BFGS",
      n_pop = n_pop,
      n_age = n_age,
      mortality_rate = log_residual,
      fitted_mortality = fitted
    )$par
    
    gamma <- c(gamma, gamma_value)
    
    log_residual <- lapply(1:n_pop, function(j) {
      log_residual[[j]] - gamma_value * fitted[[j]]
    })
    loss_vec[i] <- sum(unlist(log_residual)^2)

    
    ax[[1]][i,] <- product_fitted$ax
    bx[[1]][i,] <- product_fitted$bx
    kt[[1]][,i] <- product_fitted$kapa
    
    for(k in 1:n_pop){
      ax[[k+1]][i,] <- ratio_fitted[[k]]$ax
      bx[[k+1]][i,] <- ratio_fitted[[k]]$bx
      kt[[k+1]][,i] <- ratio_fitted[[k]]$kapa
    }
    
    # new stopping criteria
    if (i > 1) {
      if (abs(loss_vec[i] - loss_vec[i-1]) < epsilon) {
        iteration <- i
        break
      }
    }
  } # end loop
  
  
  for(z in 1:iteration){
    for(g in 1:n_pop){
      for(h in 1:n_age){
        fitted_value[[g]][,h] <- fitted_value[[g]][,h] +
          gamma[z] * (ax[[1]][z,h] + bx[[1]][z,h]*kt[[1]][,z] +
                        ax[[g+1]][z,h] + bx[[g+1]][z,h]*kt[[g+1]][,z])
      }
    }
  }
  
  fitted_value <- lapply(fitted_value, exp)
  return(list(
    ax = ax,
    kapa = kt,
    bx = bx,
    iteration = iteration,
    gamma = gamma,
    fitted_value = fitted_value,
    loss = loss_vec[1:iteration]
  ))
}






6. Forecast performance based on MASE
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
