# The prediction interval based on simulation
simulation_pi <- function(n_pop, n_age, horizon, fitted_output, train_regressor, test_regressor,
                  flip_indicator, simulation_iter = 100, gamma, quantile = c(0.025, 0.975), seed) {
  set.seed(seed)
  ax <- fitted_output$ax
  bx <- fitted_output$bx
  kt <- fitted_output$kapa
  
  n_kt <- ncol(kt[[1]])
  
  # Initialize kt forecast: list of (pop+1), each is array [horizon, n_kt, simulation_iter]
  kt_forecast <- vector("list", n_pop + 1)
  for(d in 1:(n_pop + 1)) {
    kt_forecast[[d]] <- array(0, dim = c(horizon, n_kt, simulation_iter))
  }
  
  # Simulate each kt once per simulation
  for(d in 1:(n_pop + 1)) {
    for(z in 1:n_kt) {       
      kt_series <- kt[[d]][, z]
      for(f in 1:simulation_iter) {
        fit.gbll.kt <- auto.arima(kt_series, xreg = train_regressor)
        
        # Decide whether to use xreg
        if(length(coef(fit.gbll.kt)) < 4) {
          fit.gbll.kt <- auto.arima(kt_series)
          kt_forecast[[d]][, z, f] <- simulate(fit.gbll.kt, nsim = horizon, future = TRUE)
        } else if(coeftest(fit.gbll.kt)["xreg1", "Pr(>|z|)"] < 0.05 || 
                  coeftest(fit.gbll.kt)["xreg2", "Pr(>|z|)"] < 0.05 ||
                  coeftest(fit.gbll.kt)["xreg3", "Pr(>|z|)"] < 0.05 ||
                  coeftest(fit.gbll.kt)["xreg4", "Pr(>|z|)"] < 0.05) {
          kt_forecast[[d]][, z, f] <- simulate(fit.gbll.kt, nsim = horizon, future = TRUE, xreg = test_regressor)
        } else {
          fit.gbll.kt <- auto.arima(kt_series)
          kt_forecast[[d]][, z, f] <- simulate(fit.gbll.kt, nsim = horizon, future = TRUE)
        }
      }
    }
  }
  
  # Initialize forecasted mortality log
  forecasted_mortality_log <- vector("list", n_pop)
  for(g in 1:n_pop) {
    forecasted_mortality_log[[g]] <- array(0, dim = c(horizon, n_age, simulation_iter))
  }
  
  # Construct forecasted mortality for each simulation
  for(sim in 1:simulation_iter) {
    for(g in 1:n_pop) {
      for(h in 1:n_age) {
        forecasted_mortality_log[[g]][, h, sim] <- 0
        for(z in 1:n_kt) {
          forecasted_mortality_log[[g]][, h, sim] <- forecasted_mortality_log[[g]][, h, sim] +
            gamma[z] * (
              ax[[1]][z, h] + bx[[1]][z, h] * kt_forecast[[1]][, z, sim] +   # common kt
              ax[[g + 1]][z, h] + bx[[g + 1]][z, h] * kt_forecast[[g + 1]][, z, sim]  # country-specific
            )
        }
      }
    }
  }
  
  # Convert to actual mortality
  forecasted_mortality <- forecasted_mortality_log
  for(g in 1:n_pop) {
    if(flip_indicator[g] == 0) {
      forecasted_mortality[[g]] <- exp(forecasted_mortality[[g]])
    } else {
      forecasted_mortality[[g]] <- 1 / exp(forecasted_mortality[[g]])
    }
  }
  
  # Compute prediction intervals
  lower <- vector("list", n_pop)
  upper <- vector("list", n_pop)
  for(g in 1:n_pop) {
    lower[[g]] <- apply(forecasted_mortality[[g]], c(1, 2), function(x) quantile(x, quantile[1]))
    upper[[g]] <- apply(forecasted_mortality[[g]], c(1, 2), function(x) quantile(x, quantile[2]))
  }
  
  return(list(lower = lower, upper = upper))
}
