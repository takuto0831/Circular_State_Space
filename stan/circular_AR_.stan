data{
  // data paramter
  int N; // sample size
  int P; // AR(P) 
  vector[N] theta; // data
}

parameters{
  real alpha_0;
  row_vector[P] alpha_1;
  real<lower=0> sigma;
}

model{
  // パラメータの事前分布の分散の値も入力データとする
  alpha_0 ~ normal(0,100); sigma ~ normal(0,100);
  for(i in 1:P){
    alpha_1[i] ~ normal(0,100); // N(0,100)
  }
  for(n in (1+P):N){
    vector[P] pre_theta;
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    target += normal_lpdf(theta[n]|alpha_0 + ( alpha_1 * pre_theta),sigma);
  }
}

generated quantities{
  vector[N-P] log_likelihood;
  vector[N-P] theta_new;
  for(n in (1+P):N){
    vector[P] pre_theta;
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    log_likelihood[n-P] = normal_lpdf(theta[n]|alpha_0 + ( alpha_1 * pre_theta),sigma);
    theta_new[n-P] = normal_rng(alpha_0 + ( alpha_1 * pre_theta),sigma);
  } 
}
