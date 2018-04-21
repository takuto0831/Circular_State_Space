CircularRegSim <- function(theta_vec, lag, alpha, Beta){
  #theta_vec:Angle data, alpha:(2,1), Beta:(2,2,lag)
  len <- length(theta_vec) 
  d <- matrix(nrow = 2,ncol = len - lag) 
  
  for(i in 1:len-lag){
    tmp <- matrix(c(0,0),nrow = 2)
    for(j in 1:lag){
      tmp = tmp + beta[,,j] %*% matrix(c(cos(theta_vec[i-j]),sin(theta_vec[i-j])),nrow = 2)
    }
    d[,i] = alpha + tmp 
  }
}
