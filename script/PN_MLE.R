#################### Projected Normal AutoRegressive Process (p=1) ################

# Projected Normal AutoRegressive Process
l <- function(arg){
  # data = theta, arg = c(alpha_0, alpha_1, Sigma)
  arg <- matrix(arg,nrow = 2) # データを行列化
  likelihood <- c(0)
  for(i in 2:10){
    u = matrix(c(cos(data[i]),sin(data[i])),ncol=1) 
    mu = arg[,1] + arg[,2:3] %*% matrix(c(cos(data[i-1]),sin(data[i-1])),ncol=1)
    A = t(u) %*% solve(arg[,4:5]) %*% u
    B = t(u) %*% solve(arg[,4:5]) %*% mu
    C = (-1/2) * (t(mu) %*% solve(arg[,4:5]) %*% mu)
    tmp = B/sqrt(A)
    # 尤度を計算する (-1を乗じることで, 最大化問題を最小化問題にする)
    likelihood <- append(likelihood,
                         # -log(A) - 0.5*log(det(arg[,4:5])) + C + log(1+(tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1)))
                         log(A) + 0.5*log(det(arg[,4:5])) - C - log(1+(tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1)))
    )
  }
  sum(likelihood)
}

ll <- function(arg){
  # data = theta, arg = c(alpha_0, alpha_1, tau,rho)
  arg <- matrix(arg,nrow = 2) # データを行列化
  Sigma <- matrix(c(arg[1,4]^2, arg[1,4]*arg[2,4], 
                    arg[1,4]*arg[2,4],1),nrow = 2)
  likelihood <- c(0)
  for(i in 2:2){
    u = matrix(c(cos(data[i]),sin(data[i])),ncol=1) 
    mu = arg[,1] + arg[,2:3] %*% matrix(c(cos(data[i-1]),sin(data[i-1])),ncol=1)
    A = t(u) %*% solve(Sigma) %*% u
    B = t(u) %*% solve(Sigma) %*% mu
    C = (-1/2) * (t(mu) %*% solve(Sigma) %*% mu)
    tmp = B/sqrt(A)
    # 尤度を計算する (-1を乗じることで, 最大化問題を最小化問題にする)
    likelihood <- append(likelihood,
                         # -log(A) - 0.5*log(det(Sigma)) + C + log(1+(tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1)))
                         log(A) + 0.5*log(det(Sigma)) - C - log(1+(tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1)))
    )
    browser()
  }
  #browser()
  sum(likelihood)
}

# 不等式制約 for l
inequalityConstraint <- function(arg){
  arg <- matrix(arg,nrow = 2) # データを行列化
  tmp1 <- eigen(arg[,4:5])$values
  tmp2 <- det(arg[,4:5]) %>% abs()
  c(tmp1,tmp2) # 逆行列を持つ
}

ineq.lower <- c(0.00001,0.00001,0.00001)
ineq.upper <- c(1000,1000,1000)

# 不等式制約 for ll
inequalityConstraint <- function(arg){ arg }; 
ineq.lower <-c(-100,-100,-100,-100,-100,-100,0.001,-1); 
ineq.upper <- c(100,100,100,100,100,100,100,1)

# 結果のベクトルを抽出する
ans <- function(solution){
  tmp <- solution$pars %>% matrix(nrow=2)
  Ans <- c()
  for(i in 2:num){
    ans <- tmp[,1] + (tmp[,2:3] %*% matrix(c(cos(data[i-1]),sin(data[i-1])),ncol=1)) + mvrnorm(1,c(0,0),tmp[,4:5])
    Ans <- rbind(Ans,t(ans))
  }
  return(Ans)
}

#######################  main 関数  ##########################
library(dplyr)
library(tidyverse)
library(circular)

# dataset we report the wind direction recorded every day from January 29, 2001 to March 31, 2001 from 3.00am to 4.00am included.
data(wind) 
wind_data <- wind %>%  
  data.frame(t = seq(1,310,1), theta_real=.) %>%  # put label
  mutate(cos_real = cos(theta_real), sin_real=sin(theta_real))  

#################### パラメータ初期値 #####################
alpha_0 <- matrix(c(0,0),nrow=2) # const
alpha_1 <- matrix(c(0.3,0,0,0.2),nrow=2) #mu_0
# Sigma <- matrix(c(0.1,0,0.1,0.1),nrow=2) #sigma
tau <- 0.01; rho <- 0
data <- wind_data$theta_real # analysis
#data <- Sim_data$theta # simulation
num <- length(data)
# 推定したいパラメータをまとめる
# arg <- c(alpha_0,alpha_1,Sigma)
arg <- c(alpha_0,alpha_1,tau,rho)

############## Rsolnp を用いて最適化 ######################
library(Rsolnp)
solution <- solnp(arg, fun = ll)
solution <- solnp(arg, fun = ll, ineqfun = inequalityConstraint, ineqLB =ineq.lower,ineqUB = ineq.upper) #制約付き問題data_ans <- ans(solution) # パラメータ推定値を用いて, 予測

############### 結果の出力 ################
data_ans %>% 
  as_data_frame() %>% 
  mutate(t = seq(1,length(data_ans)/2,1)) %>% 
  tidyr::gather(key = "variable",value  = value, V1,V2) %>% 
  ggplot(aes(x=t,y=value,colour=variable)) +
  geom_line() + 
  theme_economist() +
  labs(title = "time series circular plot")
