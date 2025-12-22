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


