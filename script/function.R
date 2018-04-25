# シミュレーションデータを生成する関数
CircularRegSim <- function(theta_start_vec, len, alpha, beta, Eps){
  #theta_vec:Angle data, alpha:(2,1), beta:(2,2,lag)
  lag <- length(theta_start_vec)
  mu <- matrix(c(0,0),nrow = 2)
  d <- matrix(nrow = 3, ncol = len) 
  
  # 初期値を行列dに格納する.
  for(i in 1:lag){
    d[1,i] = i # 次数
    d[2,i] = cos(theta_start_vec[i])
    d[3,i] = sin(theta_start_vec[i])
  }
  for (i in seq(lag+1,len,1)){
    tmp <- matrix(c(0,0),nrow = 2)
    for(j in 1:lag){
      tmp = tmp + (beta[,,j] %*% d[2:3,i-lag])
    }
    d[1,i] = i
    d[2:3,i] = alpha + tmp + as.matrix(mvrnorm(1,mu,Eps))
  }
  return(d)
}

# 正規乱数のパラメータより, 分散共分散行列を計算する関数
GeneratingEps <- function(sigma_c, sigma_s, rho){
  matrix(c(sigma_c^2,rho*sigma_c*sigma_s,rho*sigma_c*sigma_s,sigma_s^2),nrow=2) %>% 
    return()
}