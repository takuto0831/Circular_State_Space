functions{
  real circular_ar_lpdf(real theta, int P, vector pre_theta, real alpha_0, row_vector alpha_1, real sigma){
    real mu; real p;
    mu = alpha_0 + ( alpha_1 * pre_theta); // const + coef * theta 
    p = normal_lpdf(theta|mu,sigma);
    return p;
  }
}

data{
  // data paramter
  int N; // sample size
  int P; // AR(P) 
  vector<lower=-pi(),upper=pi()>[N] theta; // data
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
    target += circular_ar_lpdf(theta[n]|P,pre_theta, alpha_0,alpha_1,sigma);
  }
}

generated quantities{
  vector[N-P] log_likelihood;
  for(n in (1+P):N){
    vector[P] pre_theta;
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    log_likelihood[n-P] = circular_ar_lpdf(theta[n]|P,pre_theta, alpha_0,alpha_1,sigma);
  } 
}
