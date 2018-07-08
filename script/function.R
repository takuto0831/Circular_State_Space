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
# |Sigma| = 1 の制約を加える
GeneratingEps2 <- function(sigma_c, rho){
  sigma_s = sqrt( 1 / (sigma_c^2*( 1 - rho^2)));
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


########## Var model 用 ##########
ans_stan_p <- function(fit,P,data){
  # 得られたモデルを抽出する
  fit_ext <- rstan::extract(fit,permuted=T); 
  # Const parameter
  alpha_0 = matrix(c(fit_ext$alpha_0[,1] %>% mean(),fit_ext$alpha_0[,2] %>% mean()),ncol=1);
  # alpha_0_low = matrix(c(fit_ext$alpha_0[,1] %>% quantile(0.025),fit_ext$alpha_0[,2] %>% quantile(0.025)),ncol=1);
  # alpha_0_up = matrix(c(fit_ext$alpha_0[,1] %>% quantile(0.975),fit_ext$alpha_0[,2] %>% quantile(0.975)),ncol=1);
  
  # Beta parameter
  alpha_1 <- alpha_1_low <- alpha_1_up <- c();
  for(i in 1:P){
    tmp <- matrix(c(fit_ext$alpha_1[,1,(2*i-1)] %>% mean(),fit_ext$alpha_1[,1,2*i] %>% mean(),
                    fit_ext$alpha_1[,2,(2*i-1)] %>% mean(),fit_ext$alpha_1[,2,2*i] %>% mean()),byrow=T,2,2)
    alpha_1 <- cbind(alpha_1, tmp)
    # tmp_low <- matrix(c(fit_ext$alpha_1[,1,(2*i-1)] %>% quantile(0.025),fit_ext$alpha_1[,1,2*i] %>% quantile(0.025),
    #                     fit_ext$alpha_1[,2,(2*i-1)] %>% quantile(0.025),fit_ext$alpha_1[,2,2*i] %>% quantile(0.025)),byrow=T,2,2)
    # alpha_1_low <- cbind(alpha_1_low, tmp_low)
    # tmp_up <- matrix(c(fit_ext$alpha_1[,1,(2*i-1)] %>% quantile(0.975),fit_ext$alpha_1[,1,2*i] %>% quantile(0.975),
    #                    fit_ext$alpha_1[,2,(2*i-1)] %>% quantile(0.975),fit_ext$alpha_1[,2,2*i] %>% quantile(0.975)),byrow=T,2,2)
    # alpha_1_up <- cbind(alpha_1_up, tmp_up)
  }
  # Variance-Covariance matrix
  Sigma = matrix(c(fit_ext$sigma[,1,1] %>% mean(),fit_ext$sigma[,1,2] %>% mean(),
                   fit_ext$sigma[,2,1] %>% mean(),fit_ext$sigma[,2,2] %>% mean()),byrow=T,2,2)
  # Sigma_low = matrix(c(fit_ext$sigma[,1,1] %>% quantile(0.025),fit_ext$sigma[,1,2] %>% quantile(0.025),
  #                      fit_ext$sigma[,2,1] %>% quantile(0.025),fit_ext$sigma[,2,2] %>% quantile(0.025)),byrow=T,2,2)
  # Sigma_up = matrix(c(fit_ext$sigma[,1,1] %>% quantile(0.975),fit_ext$sigma[,1,2] %>% quantile(0.975),
  #                     fit_ext$sigma[,2,1] %>% quantile(0.975),fit_ext$sigma[,2,2] %>% quantile(0.975)),byrow=T,2,2)
  
  # Estimate condition mean 
  num = length(data) 
  mu_hat <- mu_hat_low <- mu_hat_up <- matrix(0, ncol=2, nrow=(num-P) )
  for(i in (1+P):num){
    pre <- c(); # p期前のcos(theta), sin(theta)を格納する
    for (k in 1:P) pre <- rbind(pre,matrix(c(cos(data[i-k]),sin(data[i-k])),ncol=1));
    mu_hat[i-P,] <- alpha_0 + ( alpha_1 %*% pre )
    #mu_hat_low[i-P,] <- alpha_0_low + ( alpha_1_low %*% pre )
    #mu_hat_up[i-P,] <- alpha_0_up + ( alpha_1_up %*% pre )
  }
  # Estimate moment function
  fn.sin<- function(x,mu,Sigma) sin(x)*v.dPnCircular_dens(theta=x, mu, Sigma) # sin func
  fn.cos<- function(x,mu,Sigma) cos(x)*v.dPnCircular_dens(theta=x, mu, Sigma) # cos func
  fn.theta <- function(x,mu,Sigma) x*v.dPnCircular_dens(theta=x, mu, Sigma) # theta func
  
  ########################## Estimate theta value through sin and cos value ################################
  # Estimate sin and cos value
  len <- dim(mu_hat)[1];
  sin.mom <- cos.mom <- sin.mom_low <- cos.mom_low <- sin.mom_up <- cos.mom_up <- c();
  for(i in 1:len){
    cos.mom[i] <- integrate(fn.cos, lower=-pi, upper=pi, mu=mu_hat[i,], Sigma=Sigma)$value
    sin.mom[i] <- integrate(fn.sin, lower=-pi, upper=pi, mu=mu_hat[i,], Sigma=Sigma)$value
    # sin.mom_low[i] <- integrate(fn.sin, lower=-pi, upper=pi, mu=mu_hat_low[i,], Sigma=Sigma_low)$value
    # cos.mom_low[i] <- integrate(fn.cos, lower=-pi, upper=pi, mu=mu_hat_low[i,], Sigma=Sigma_low)$value
    # sin.mom_up[i] <- integrate(fn.sin, lower=-pi, upper=pi, mu=mu_hat_up[i,], Sigma=Sigma_up)$value
    # cos.mom_up[i] <- integrate(fn.cos, lower=-pi, upper=pi, mu=mu_hat_up[i,], Sigma=Sigma_up)$value
  }
  # Estimate predict theta value
  pred <- atan2(sin.mom,cos.mom)
  # pred_low <- atan2(sin.mom_low,cos.mom_low)
  # pred_up <- atan2(sin.mom_up,cos.mom_up)
  
  ########################## Estimate theta value directly ################################
  # len <- dim(mu_hat)[1]; 
  # pred <- pred_low <- pred_up <- c();
  # for(i in 1:len){
  #  pred[i] <- integrate(fn.theta, lower=-pi, upper=pi, mu=mu_hat[i,], Sigma=Sigma)$value
  #  pred_low[i] <- integrate(fn.theta, lower=-pi, upper=pi, mu=mu_hat_low[i,], Sigma=Sigma_low)$value
  #  pred_up[i] <- integrate(fn.theta, lower=-pi, upper=pi, mu=mu_hat_up[i,], Sigma=Sigma_up)$value
  # }
  
  # return(data.frame(predict= c(rep(NA,P),pred),lower = c(rep(NA,P),pred_low),upper = c(rep(NA,P),pred_up)))  
  return(data.frame(predict= c(rep(NA,P),pred)))
}

########## ar model 用 ##########
ans_stan_p2 <- function(fit,P,data){
  # 得られたモデルを抽出する
  fit_ext <- rstan::extract(fit,permuted=T); 
  # Const parameter
  alpha_0 = fit_ext$alpha_0 %>% mean();
  # Beta parameter
  alpha_1 <- c();
  for(i in 1:P){
    tmp <- fit_ext$alpha_1[,i] %>% mean(); alpha_1 <- cbind(alpha_1, tmp); 
  }
  # Variance-Covariance matrix
  Sigma = fit_ext$sigma %>% mean()
  
  # Estimate condition mean 
  num = length(data); mu_hat <- c() 
  for(i in (1+P):num){
    tmp <- alpha_0 + ( alpha_1 %*% data[(i-1):(i-P)] )
    mu_hat <- cbind(mu_hat, tmp)
  }
  # Estimate moment function
  fn.theta <- function(k,mu,Sigma) k*dnorm(x=k, mu, Sigma) # theta func

  ########################## Estimate theta value through sin and cos value ################################
  # Estimate sin and cos value
  len <- length(mu_hat); theta.mom <- c();
  for(i in 1:len){
    theta.mom[i] <- integrate(fn.theta, lower=-Inf, upper=Inf, mu=mu_hat[i], Sigma=Sigma)$value
  }
  # Estimate predict theta value
  return(data.frame(predict= c(rep(NA,P),theta.mom)))
}
# output predict value for VAR
pred_value <- function(fit,p,dat,who=2){
  # who:1 -> AR model, who:2 -> VAR model
  if(who == 1){pred <- ans_stan_p2(fit, P=p, data=dat)} # ar model
  else{pred <- ans_stan_p(fit, P=p, data=dat)} # var model
  # model のRMSE 算出するコード
  data.frame(real=dat,pred) %>% 
    mutate(dif = predict - real) %>% 
    mutate(dif_ = dplyr::if_else(dif < -pi, dif + 2*pi,
                                 dplyr::if_else(dif > pi, dif - 2*pi, dif))) %>% # 誤差は必ず -pi ~ pi に含まれる
    mutate(eps = dif_^2) %>% 
    summarise(ans = sqrt(mean(eps,na.rm = TRUE))) %>% 
    {paste("RMSE of", fit@model_name, "model.",sprintf("lag = %d, RMSE  = %f",p,.),sep = " ")} %>% 
    print() # Root Mean Square Error
  # 予測値, 真値, 信用区間を表示する
  data.frame(real=dat,pred) %>% 
    mutate(index = row_number()) %>%
    ggplot(aes(x=index)) +
    # geom_ribbon(aes(ymin=lower,ymax=upper),fill="skyblue",alpha=0.6) +
    geom_line(aes(y=real, colour= "real")) + 
    geom_line(aes(y=predict, colour= "predict")) +
    facet_zoom(x = index >= (length(dat) -50)) + 
    scale_colour_manual("",
                        values = c("real"="black","predict"="blue")
    ) +
    labs(y=expression(theta))
}


# output PN distribution in arbitrary index
PN_dist_pred <- function(fit,index_vec,sample_num,data,P){
  fit_ext <- rstan::extract(fit,permuted=T);
  # Const parameter
  alpha_0 <- fit_ext$alpha_0 %>% t()
  # Beta parameter
  alpha_1 <- fit_ext$alpha_1 # 160000, 2, 2*p
  # Variance-Covariance matrix
  Sigma <- fit_ext$sigma # 16000,2,2
  # 空配列の定義
  ans_box <- data.frame(id = seq(1,sample_num))
  # 積分関数定義
  fn.sin<- function(x,mu,Sigma) sin(x)*v.dPnCircular_dens(theta=x, mu, Sigma) # sin func
  fn.cos<- function(x,mu,Sigma) cos(x)*v.dPnCircular_dens(theta=x, mu, Sigma) # cos func
  # Estimate condition mean
  for (i in 1:length(index_vec)) {
    pre <- c(); # p期前のcos(theta), sin(theta)を格納する
    cos.mom <- sin.mom <- c();
    for (k in 1:P) pre <- rbind(pre,matrix(c(cos(data[index_vec[i] - k]),sin(data[index_vec[i] - k])),ncol=1));
    for (k in 1:sample_num) {
      mu_hat <- alpha_0[,k] + ( alpha_1[k,,] %*% pre )
      cos.mom[k] <- integrate(fn.cos, lower=-pi, upper=pi, mu=mu_hat, Sigma=Sigma[k,,])$value
      sin.mom[k] <- integrate(fn.sin, lower=-pi, upper=pi, mu=mu_hat, Sigma=Sigma[k,,])$value
    }
    pred <- atan2(sin.mom,cos.mom)
    ans_box <- cbind(ans_box,pred)
    colnames(ans_box)[i+1] <- paste("index",index_vec[i],sep = "")
  }
  return(ans_box)
}

DIC_func_AR <- function(fit,data,dev,P)
{
  ## mean.Deviance
  # data = theta, dev = log PN dist
  mean.Deviance <- extract_log_lik(fit,"log_likelihood") %>% # nrow = 16000, ncol = n - P,の対数尤度 
    apply(1, sum) %>% # 各サンプルでの尤度の総和
    mean() # 各サンプルでの尤度の総和の平均
  mean.Deviance <- -2 * mean.Deviance # -2倍
  
  # 得られたモデルを抽出する
  fit_ext <- rstan::extract(fit,permuted=T); 
  # Const parameter
  alpha_0 = fit_ext$alpha_0 %>% mean();
  # Beta parameter
  alpha_1 <- c();
  for(i in 1:P){
    tmp <- fit_ext$alpha_1[,i] %>% mean(); alpha_1 <- cbind(alpha_1, tmp); 
  }
  # Variance-Covariance matrix
  Sigma = fit_ext$sigma %>% mean()
  
  # Estimate condition mean 
  num = length(data); mu_hat <- c() 
  for(i in (1+P):num){
    tmp <- alpha_0 + ( alpha_1 %*% data[(i-1):(i-P)] )
    mu_hat <- cbind(mu_hat, tmp)
  }
  len <- length(mu_hat); Deviance.mean <- 0; # パラメータの平均値を用いて, 算出した尤度の総和
  for(i in 1:len){
    Deviance.mean <- Deviance.mean + (-2 * as.numeric(dev(data[i+P], mu_hat[i], Sigma,log = T)))
  }
  pD <- mean.Deviance - Deviance.mean
  return(c(P = P,mean.Dev=mean.Deviance, Dev.mean=Deviance.mean, pD=mean.Deviance-Deviance.mean, DIC=pD+mean.Deviance) )
}

DIC_func_VAR <- function(fit,data,dev,P)
{
  ## mean.Deviance
  # data = theta, dev = log PN dist
  mean.Deviance <- extract_log_lik(fit,"log_likelihood") %>% # nrow = 16000, ncol = n - P,の対数尤度 
    apply(1, sum) %>% # 各サンプルでの尤度の総和
    mean() # 各サンプルでの尤度の総和の平均
  mean.Deviance <- -2 * mean.Deviance # -2倍
  
  ## Deviance.mean
  fit_ext <- rstan::extract(fit,permuted=T); 
  # Const parameter
  alpha_0 = matrix(c(fit_ext$alpha_0[,1] %>% mean(),fit_ext$alpha_0[,2] %>% mean()),ncol=1);
  # Beta parameter
  alpha_1 <- c();
  for(i in 1:P){
    tmp <- matrix(c(fit_ext$alpha_1[,1,(2*i-1)] %>% mean(),fit_ext$alpha_1[,1,2*i] %>% mean(),
                    fit_ext$alpha_1[,2,(2*i-1)] %>% mean(),fit_ext$alpha_1[,2,2*i] %>% mean()),byrow=T,2,2)
    alpha_1 <- cbind(alpha_1, tmp)
  }
  # Variance-Covariance matrix
  Sigma = matrix(c(fit_ext$sigma[,1,1] %>% mean(),fit_ext$sigma[,1,2] %>% mean(),
                   fit_ext$sigma[,2,1] %>% mean(),fit_ext$sigma[,2,2] %>% mean()),byrow=T,2,2)
  
  # Estimate condition mean 
  num = length(data) 
  mu_hat <- matrix(0, ncol=2, nrow=(num-P) )
  for(i in (1+P):num){
    pre <- c(); # p期前のcos(theta), sin(theta)を格納する
    for (k in 1:P) pre <- rbind(pre,matrix(c(cos(data[i-k]),sin(data[i-k])),ncol=1));
    mu_hat[i-P,] <- alpha_0 + ( alpha_1 %*% pre )
  }
  len <- dim(mu_hat)[1]; Deviance.mean <- 0; # パラメータの平均値を用いて, 算出した尤度の総和
  for(i in 1:len){
    Deviance.mean <- Deviance.mean + (-2 * as.numeric(dev(data[i+P], mu_hat[i,], Sigma)))
  }
  pD <- mean.Deviance - Deviance.mean
  return(c(P = P,mean.Dev=mean.Deviance, Dev.mean=Deviance.mean, pD=mean.Deviance-Deviance.mean, DIC=pD+mean.Deviance) )
}

# 可視化用
stan_ac_label <- function (object, pars, label_set, include = TRUE, unconstrain = FALSE, 
                           inc_warmup = FALSE, nrow = NULL, ncol = NULL, ..., separate_chains = FALSE, 
                           lags = 25, partial = FALSE) 
{
  plot_data <- .make_plot_data(object, pars, include, inc_warmup, 
                               unconstrain)
  # dots <- .add_aesthetics(list(...), c("size", "color", "fill"))
  thm <- .rstanvis_defaults$theme
  dat_args <- list(dat = plot_data$samp, lags = lags, partial = partial)
  dat_fn <- ifelse(plot_data$nparams == 1, ".ac_plot_data", ".ac_plot_data_multi")
  ac_dat <- do.call(dat_fn, dat_args)
  # main 
  # dots$position <- "dodge"; dots$stat <- "summary"; dots$fun.y <- "mean"
  y_lab <- paste("Avg.", if (partial) "partial", "autocorrelation")
  ac_labs <- labs(x = "Lag", y = y_lab)
  y_scale <- scale_y_continuous(labels = seq(0, 1, 0.25))
  base <- ggplot(ac_dat, aes_string(x = "lag", y = "ac"))
  graph <- base + 
    geom_bar(stat = "kentity", fill = "blue") +
    y_scale + ac_labs + thm +
    facet_wrap(~parameters, nrow = nrow, ncol = ncol, scales = "free_x",labeller = label_set) +
    theme(strip.text.x = element_text(size = 12))
  return(graph)
}