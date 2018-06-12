# シミュレーションデータを生成する関数
CircularRegSim <- function(theta_start_vec, len, alpha, beta, Eps){
  #theta_vec:Angle data, alpha:(2,1), beta:(2,2,lag)
  lag <- length(theta_start_vec)
  mu <- matrix(c(0,0),nrow = 2)
  d <- matrix(nrow = 3, ncol = len) 
  set.seed(31) # 乱数シード
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
    eps <- mvrnorm(1,mu,Eps)
    d[2:3,i] = alpha + tmp + eps/c(sqrt(t(eps) %*% eps ))
  }
  return(d)
}

# 予測データを生成する関数
CircularRegPred <- function(vec, lag, alpha, beta, Eps){
  #theta_vec:Angle data, alpha:(2,1), beta:(2,2,lag)
  len <- length(vec)
  mu <- matrix(c(0,0),nrow = 2)
  d <- matrix(nrow = 3, ncol = len) 
  
  # 初期値を行列dに格納する.
  for(i in 1:lag){
    d[1,i] = i # 次数
    d[2,i] = cos(vec[i])
    d[3,i] = sin(vec[i])
  }
  for (i in seq(lag+1,len,1)){
    tmp <- matrix(c(0,0),nrow = 2)
    for(j in 1:lag){
      tmp = tmp + (beta[[j]] %*% d[2:3,i-lag])
    }
    d[1,i] = i
    d[2:3,i] = alpha + tmp + as.matrix(mvrnorm(1,mu,Eps))
  }
  return(d %>% 
           t() %>% 
           as.data.frame() %>%
           rename(p=V1, cos=V2, sin=V3))
}

# パラメータを用いた予測
CircularRegPred_parameter <- function(theta,tmp,lag){
  alpha <- Bcoef(tmp) %>% as_data_frame() %>% {as.matrix(.$const)}
  beta <- Acoef(tmp) 
  Eps <- summary(tmp)$covres %>% matrix(nrow = 2)
  CircularRegPred(theta,lag,alpha,beta,Eps) %>% 
    return()
} 

# # Projected normal distribution density funciton
# PnCircular_dens <- function(theta,mu,Sigma){
#   u = matrix(c(cos(theta),sin(theta)),ncol=1)
#   A = t(u) %*% solve(Sigma) %*% u
#   B = t(u) %*% solve(Sigma) %*% mu
#   C = (-1/2) * (t(mu) %*% solve(Sigma) %*% mu)
#   tmp = B/sqrt(A)
#   p = (1/(2*pi*A*sqrt(det(Sigma)))) * exp(C) * 
#     (1 + tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1))
#   return(p)
# }

# Projected Normal distribution log likelihood 
PnCircular_log <- function(theta, mu, Sigma){
  u = matrix(c(cos(theta),sin(theta)),ncol=1)
  A = t(u) %*% solve(Sigma) %*% u
  B = t(u) %*% solve(Sigma) %*% mu
  C = (-1/2) * (t(mu) %*% solve(Sigma) %*% mu)
  tmp = B/sqrt(A)
  p = -log(A) -0.5*log(det(Sigma)) + C + log(1+(tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1)));  
  return(p)
}

# 正規乱数のパラメータより, 分散共分散行列を計算する関数
GeneratingEps <- function(sigma_c, sigma_s, rho){
  matrix(c(sigma_c^2,rho*sigma_c*sigma_s,rho*sigma_c*sigma_s,sigma_s^2),nrow=2) %>% 
    return()
}

# Projected Normal Circular Plobability Density 
dPnCircular_dens <- function(theta,mu,Sigma){
  u = matrix(c(cos(theta),sin(theta)),ncol=1)
  A = t(u) %*% solve(Sigma) %*% u
  B = t(u) %*% solve(Sigma) %*% mu
  C = (-1/2) * (t(mu) %*% solve(Sigma) %*% mu)
  tmp = B/sqrt(A)
  p = as.numeric((1/(2*pi*A*sqrt(det(Sigma)))) * exp(C) *
                   (1 + tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1)))
  return(p)
}
# Vectorize 関数
v.dPnCircular_dens <- Vectorize(dPnCircular_dens, "theta") 

# 結果のベクトルを抽出 VAR(p)
ans_stan_p <- function(fit,P,num,data){
  # 得られたモデルを抽出する
  fit_ext <- rstan::extract(fit,permuted=T)
  # 格納する変数を用意する
  alpha_1 <- c()
  # Const parameter
  alpha_0 = matrix(c(fit_ext$alpha_0[,1] %>% mean(),fit_ext$alpha_0[,2] %>% mean()),ncol=1)
  # Beta parameter
  for(i in 1:P){
    tmp <- matrix(c(fit_ext$alpha_1[,1,(2*i-1)] %>% mean(),fit_ext$alpha_1[,1,2*i] %>% mean(),
                    fit_ext$alpha_1[,2,(2*i-1)] %>% mean(),fit_ext$alpha_1[,2,2*i] %>% mean()),byrow=T,2,2)
    alpha_1 <- cbind(alpha_1, tmp)
  }
  # Variance-Covariance matrix
  Sigma_hat = matrix(c(fit_ext$sigma[,1,1] %>% mean(),fit_ext$sigma[,1,2] %>% mean(),
                       fit_ext$sigma[,2,1] %>% mean(),fit_ext$sigma[,2,2] %>% mean()),byrow=T,2,2)
  # Sigma_hat = matrix(c(1,0,0,1),ncol=2) # for model3
  # Estimate condition mean 
  mu_hat <- matrix(0, ncol=2, nrow=(num-P) )
  for(i in (1+P):num){
    pre <- c(); # p期前のcos(theta), sin(theta)を格納する
    for (k in 1:P) {
      tmp <- matrix(c(cos(data[i-k]),sin(data[i-k])),ncol=1)
      pre <- rbind(pre,tmp)
    }
    mu_hat[i-P,] <- alpha_0 + ( alpha_1 %*% pre )
  }
  # Estimate moment
  fn.sin<-function(x,mu,Sigma) sin(x)*v.dPnCircular_dens(theta=x, mu, Sigma) # sin func
  fn.cos<-function(x,mu,Sigma) cos(x)*v.dPnCircular_dens(theta=x, mu, Sigma) # cos func
  # Estimate parameter
  len <- dim(mu_hat)[1]; sin.mom <- cos.mom<- c();
  # 
  for(i in 1:len){
    sin.mom[i] <- integrate(fn.sin, lower=-pi, upper=pi,
                            mu=mu_hat[i,], Sigma=Sigma_hat)$value
    cos.mom[i] <- integrate(fn.cos, lower=-pi, upper=pi,
                            mu=mu_hat[i,], Sigma=Sigma_hat)$value
  }
  pred <- atan2(sin.mom,cos.mom)
  return(pred)  
}
# output predict value
pred_value <- function(fit,p,dat){
  pred <- ans_stan_p(fit, P=p, num = length(dat),data=dat) # パラメータ推定値を用いて, 予測
  matplot(cbind(dat[-c(1:p)],pred),type="l")
}
