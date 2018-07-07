functions{
  // state model
  real circular_state_lpdf(vector u, int P, vector pre_theta, vector alpha_0, matrix alpha_1, matrix sigma){
    vector[2*P] tmp; vector[2] mu; real p;
    for(k in 1:P){
      tmp[2*k-1] = cos(pre_theta[k]); 
      tmp[2*k] = sin(pre_theta[k]);
    }
    mu = alpha_0 + ( alpha_1 * tmp); 
    p = multi_normal_lpdf(u|mu,sigma);
    return p;
  }
  // observation model
  real circular_obs_lpdf(real theta, int P, vector pre_theta, vector alpha_0, matrix alpha_1, matrix sigma){
    vector[2] u; vector[2] mu; vector[2*P] tmp; 
    real A; real B; real C; real D; real p;
    for(k in 1:P){
      tmp[2*k-1] = cos(pre_theta[k]); 
      tmp[2*k] = sin(pre_theta[k]);
    }
    mu = alpha_0 + ( alpha_1 * tmp); 
    u[1] = cos(theta); u[2] = sin(theta);
    A = quad_form(inverse_spd(sigma), u); B = u' * inverse_spd(sigma) * mu;
    C = (-0.5) * quad_form(inverse_spd(sigma), mu); D = B/sqrt(A);
    p = - log(A) - 0.5*log(determinant(sigma)) + C + log(1+(D * (normal_cdf(D,0,1)/ (exp(-D^2 /2)/sqrt(2*pi()))))); 
    return p;
  }
}

data{
  // data paramter
  int N; // sample size
  int P; // VAR(P) 
  vector<lower=-pi(),upper=pi()>[N] theta; // data
}

transformed data{
  vector[2] u[N];
  for(k in 1:N){
    u[k,1] = cos(theta[k]);
    u[k,2] = sin(theta[k]);
  }
}

parameters{
  vector[2] alpha_0;
  matrix[2,2*P] alpha_1;
  cov_matrix[2] sigma;
}

model{
  // パラメータの事前分布の分散の値も入力データとする
  alpha_0[1] ~ normal(0,100); alpha_0[2] ~ normal(0,100);
  for(i in 1:2*P){
    alpha_1[1,i] ~ normal(0,100); // N(0,100)
    alpha_1[2,i] ~ normal(0,100); // N(0,100)
  }
  sigma ~ inv_wishart(4,diag_matrix(rep_vector(1,2)));
  for(n in (1+P):N){
    vector[P] pre_theta;
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    target += circular_state_lpdf(u[n]|P,pre_theta, alpha_0,alpha_1,sigma);
    target += circular_obs_lpdf(theta[n]|P,pre_theta,alpha_0,alpha_1,sigma);
  }
}

generated quantities{
  vector[N-P] log_likelihood;
  for(n in (1+P):N){
    vector[P] pre_theta;
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    log_likelihood[n-P] = circular_state_lpdf(u[n]|P,pre_theta, alpha_0,alpha_1,sigma) + 
                          circular_obs_lpdf(theta[n]|P, pre_theta, alpha_0, alpha_1, sigma);
  } 
}

