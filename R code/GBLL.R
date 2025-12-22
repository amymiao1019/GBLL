# First define the Lee-Carter model
lc <- function(data){
  ax <- apply(data,2,mean, na.rm=TRUE)
  residual_sd <- sweep(data, 2, ax)
  U <- svd(residual_sd)$u
  V <- svd(residual_sd)$v
  d <- svd(residual_sd)$d
  kapa <- d[1]*U[,1]*sum(V[,1])
  log_sd_fitted <- as.matrix(kapa)%*%t(V[,1]/sum(V[,1]))
  bx <- t(V[,1]/sum(V[,1]))
  log_fitted <- sweep(log_sd_fitted,2,-ax)
  
  output <- list(ax=ax, kapa=kapa, bx=bx, fitted=log_fitted)
  return(output)
}


# Create Fourier Regressors in the forecasting process
fourier_regressor <- function(date, train_time, horizon, K){
  # create the train and test weekly date
  date.new <- data.frame(date = date[1:train_time])
  train_weekly <-  date.new%>%
    mutate(date=as.Date(date,format='%d/%m/%Y'),
           ddate=decimal_date(date)%%1)%>%pull(ddate)
  test_weekly <- decimal_date(date[(train_time+1):(train_time+horizon)])%%1
  
  # create the fourier regressor
  X_train <- matrix(NA, train_time, 2*K)
  X_test <- matrix(NA, horizon, 2*K)
  
  for(k in 1:K){
    X_train[,(k-1)*2+1] <- sin(2*pi*k*train_weekly)
    X_train[,k*2] <- cos(2*pi*k*train_weekly)
    
    X_test[,(k-1)*2+1] <- sin(2*pi*k*test_weekly)
    X_test[,k*2] <- cos(2*pi*k*test_weekly)
  }
  
  output <- list(train=X_train, test=X_test)
  return(output)
}


# Create the objective function to tune the learning rate gamma in the GBLL framework
obj <- function(n_pop, n_age, mortality_rate, fitted_mortality, gamma){
  obj_sum <- 0
  if(n_pop == 1){
    for(j in 1:n_age){
      fitted <- fitted_mortality[,j]
      error <- mortality_rate[,j] - gamma * fitted
      sq_error <- t(error)%*%error
      obj_sum <- obj_sum + sq_error
    }
  } else{
    for(i in 1:n_pop){
      for(j in 1:n_age){
        fitted <- fitted_mortality[[i]][,j]
        error <- mortality_rate[[i]][,j] - gamma * fitted
        sq_error <- t(error)%*%error
        obj_sum <- obj_sum + sq_error
      }
    }
  }
  return(obj_sum)
}


# GBLL model
gbll_original <- function(n_pop, n_age, mortality_rate_list, train_time, horizon, max_iteration){
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
    
    gamma_value <- optim(par=1, obj, method = "BFGS", n_pop = n_pop, n_age = n_age, mortality_rate = log_residual, fitted_mortality = fitted)$par
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







# GBLL forecasting
gbll_forecasting <- function(n_pop, n_age, horizon, fitted_output, train_regressor, test_regressor, flip_indicator, test_mortality, week_indicator, h_months, iteration, gamma){
  
  # rename and store the values of ax, bx, kt
  ax <- fitted_output$ax
  bx <- fitted_output$bx
  kt <- fitted_output$kapa
  
  kt_forecast <- create_matrices_list((n_pop+1), iteration, horizon)
  
  # forecast the kt
  for(d in 1:(n_pop+1)){
    for(f in 1:iteration){
      fit.gbll.kt <- auto.arima(kt[[d]][,f], xreg = train_regressor)
      
      if(length(coef(fit.gbll.kt)) < 4){
        fit.gbll.kt <- auto.arima(kt[[d]][,f])
        kt_forecast[[d]][,f] <- forecast(fit.gbll.kt, h = horizon)$mean
      } else if(coeftest(fit.gbll.kt)["xreg1", "Pr(>|z|)"] < 0.05 || 
                coeftest(fit.gbll.kt)["xreg2", "Pr(>|z|)"] < 0.05 ||
                coeftest(fit.gbll.kt)["xreg3", "Pr(>|z|)"] < 0.05 ||
                coeftest(fit.gbll.kt)["xreg4", "Pr(>|z|)"] < 0.05){
        kt_forecast[[d]][,f] <- forecast(fit.gbll.kt, xreg = test_regressor)$mean
      } else{
        fit.gbll.kt <- auto.arima(kt[[d]][,f])
        kt_forecast[[d]][,f] <- forecast(fit.gbll.kt, h = horizon)$mean
      }
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
  forecast_error <- create_matrices_list(n_pop, 1, h_months)
  
  difference <- Map(`-`, forecasted_mortality, test_mortality)
  abs_difference <- lapply(difference, function(diff) abs(diff))
  forecast_error_individual <- Map(`/`, abs_difference, test_mortality)
  
  rowsum_forecast_error <- lapply(forecast_error_individual, function(error) rowSums(error))
  
  for(m in 1:n_pop){
    for(n in 1:h_months){
      forecast_error[[m]][n] <- sum(rowsum_forecast_error[[m]][1:week_indicator[n+1]])/(n_age*week_indicator[n+1])
    }
  }
  
  
  # return the value
  output <- list(forecasted_mortality=forecasted_mortality, error=forecast_error, individual_error = forecast_error_individual)
  return(output)
}
