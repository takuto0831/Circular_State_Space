#################### Projected Normal AutoRegressive Process (p=1) ################

# Projected Normal AutoRegressive Process (p=1)
l <- function(arg){
  # data = theta, arg = c(alpha_0, alpha_1, Sigma)
  arg = matrix(arg,nrow = 2) # データを行列化
  likelihood <- c(0)
  for(i in 2:num){
    u = matrix(c(cos(data[i]),sin(data[i])),ncol=1) 
    mu = arg[,1] + arg[,2:3] %*% matrix(c(cos(data[i-1]),sin(data[i-1])),ncol=1)
    A = t(u) %*% solve(arg[,4:5]) %*% u
    B = t(u) %*% solve(arg[,4:5]) %*% mu
    C = (-1/2) * (t(mu) %*% solve(arg[,4:5]) %*% mu)
    tmp = B/sqrt(A)
    likelihood <- append(likelihood,
                         -log(A) - 0.5*log(det(arg[,4:5])) + C + log(1+(tmp*pnorm(tmp,0,1)/dnorm(tmp,0,1)))
    )
    }
  return(-sum(likelihood))
}

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

#################### パラメータ初期値 #####################
alpha_0 <- matrix(c(0,0),nrow=2) # const
alpha_1 <- matrix(c(0.3,0,0,0.2),nrow=2) #mu_0
Sigma <- matrix(c(0.1,0,0.1,0.1),nrow=2) #sigma
#data <- wind_data$theta_real # analysis
data <- Sim_data$theta # simulation
num <- length(data)
# 推定したいパラメータをまとめる
arg <- c(alpha_0,alpha_1,Sigma)

############## Rsolnp を用いて最適化 ######################
solution <- solnp(arg, fun = l)
# solution <- solnp(arg, fun = l, eqfun = ConstFunc, eqB = eq.value ) #制約付き問題

############### 結果の出力 #################
data_ans <- ans(solution)

data_ans %>% 
  as_data_frame() %>% 
  mutate(t = seq(1,length(data_ans)/2,1)) %>% 
  tidyr::gather(key = "variable",value  = value, V1,V2) %>% 
  ggplot(aes(x=t,y=value,colour=variable)) +
  geom_line() + 
  theme_economist() +
  labs(title = "time series circular plot")


############# stan code ################
library(rstan)
d.dat<-list(N=length(data),theta=data,pre_theta = data[-1])
d.fit<-stan(file='stan/test.stan',data=d.dat,iter=500,chains=1)
