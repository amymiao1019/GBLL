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


