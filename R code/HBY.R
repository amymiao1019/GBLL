# For the HBY model, please you the data provided in the "HBY data" folder
hby_fit <- function(name, mortality_rate, exposure, order, n_pop, n_age, train_time){
  # create empty matrices list to store value
  kt <- create_matrices_list((n_pop+1), order, train_time)
  ax <- create_matrices_list((n_pop+1), n_age, order)
  bx <- create_matrices_list((n_pop+1), n_age, order)
  fitted_value <- create_matrices_list(n_pop, n_age, train_time)
  
  
  # fit HBY
  data <- read.demogdata(mortality_rate, exposure, type = "mortality", label = name, skip = 0, popskip = 0, max.mx=9000)
  HBY.fit <- coherentfdm(data, order1=order, order2=order)
  
  
  # store the value for the common trend first
  ax[[1]][1,] <- HBY.fit[["product"]][["basis"]][,1]
  for(i in 1:order){
    kt[[1]][,i] <- HBY.fit[["product"]][["coeff"]][,i+1]
    bx[[1]][i,] <- HBY.fit[["product"]][["basis"]][,i+1]
  }
  
  # store the values for the individual trend
  for(j in 1:n_pop){
    ax[[j+1]][1,] <- HBY.fit[["ratio"]][[name[j]]][["basis"]][,1]
    for(k in 1:order){
      kt[[j+1]][,k] <- HBY.fit[["ratio"]][[name[j]]][["coeff"]][,k+1]
      bx[[j+1]][k,] <- HBY.fit[["ratio"]][[name[j]]][["basis"]][,k+1]
    }
  }
  
  # fitted values
  for(z in 1:order){
    for(g in 1:n_pop){
      for(h in 1:n_age){
        fitted_value[[g]][,h] <- fitted_value[[g]][,h] + ax[[1]][z,h] + bx[[1]][z,h]*kt[[1]][,z] + ax[[g+1]][z,h] + bx[[g+1]][z,h]*kt[[g+1]][,z]
      }
    }
  }
  
  for(j in 1:n_pop){
    fitted_value[[j]] = exp(fitted_value[[j]])
  }

  
  # return the output
  output <- list(ax=ax, kapa=kt, bx=bx, fitted_value = fitted_value)
  return(output) 
}
