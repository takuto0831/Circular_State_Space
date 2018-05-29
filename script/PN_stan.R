#################### Projected Normal AutoRegressive Process (p=1) ################

# 結果のベクトルを抽出
ans_stan <- function(fit){
  fit_ext <- rstan::extract(fit,permuted=T)
  # シュミレーション値1
  alpha_0 = matrix(c(fit_ext$alpha_0[,1] %>% mean(),fit_ext$alpha_0[,2] %>% mean()),ncol=1)
  alpha_1 = matrix(c(fit_ext$alpha_1[,1,1] %>% mean(),fit_ext$alpha_1[,1,2] %>% mean(),
                     fit_ext$alpha_1[,2,1] %>% mean(),1),ncol=2)
  # mu = matrix(c(fit_ext$mu[,1] %>% mean(),fit_ext$mu[,2] %>% mean()),ncol=1)
  Sigma = matrix(c(fit_ext$sigma[,1,1] %>% mean(),fit_ext$sigma[,1,2] %>% mean(),
                   fit_ext$sigma[,2,1] %>% mean(),fit_ext$sigma[,2,2] %>% mean()),ncol=2)
  Ans <- c()
  for(i in 2:num){
    ans <- alpha_0+ (alpha_1 %*% matrix(c(cos(data[i-1]),sin(data[i-1])),ncol=1)) + mvrnorm(1,c(0,0),Sigma)
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

############# stan code ################
library(rstan)
d.dat<-list(N=length(data),theta=data)
fit<-stan(file='stan/test.stan',data=d.dat,iter=1000,chains=1) # sigma(2,2) に対して1を仮定する
fit<-stan(file='stan/test1.stan',data=d.dat,iter=4000,chains=1) # 分散共分散行列の過程
fit<-stan(file='stan/circularVAR.stan',data=d.dat,iter=1000,chains=1) # 分散共分散行列の過程
data_ans <- ans_stan(fit) # パラメータ推定値を用いて, 予測 

############### 結果の出力 ################
data_ans %>% 
  as_data_frame() %>% 
  mutate(t = seq(1,length(data_ans)/2,1)) %>% 
  tidyr::gather(key = "variable",value  = value, V1,V2) %>% 
  ggplot(aes(x=t,y=value,colour=variable)) +
  geom_line() + 
  theme_economist() +
  labs(title = "time series circular plot")

############# 結果診断 ############ 
stan_trace(fit)
