functions{
  real circular_reg_lpdf(real theta, int P, vector pre_theta, vector alpha_0, matrix alpha_1, matrix sigma){
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
    p = - log(A) - 0.5*log(determinant(sigma)) + C
    + log(1+(D * normal_cdf(D,0,1)/exp(normal_lpdf(D|0,1))));    
    return p;
  }
}
  
data{
  int N; // sample size
  int P; // VAR(P) 
  real<lower=-pi(),upper=pi()> theta[N]; // data
}

parameters{
  vector<lower=-1,upper=1>[2] alpha_0;
  matrix<lower=-1,upper=1>[2,2*P] alpha_1; // P個の係数行列
  real phi1;
  real phi2;
  real phi3;
}

transformed parameters{
  cov_matrix[2] sigma;
  sigma[1,1] = exp(phi1)^2; sigma[1,2] = tanh(phi3)*exp(phi1)*exp(phi2);
  sigma[2,1] = tanh(phi3)*exp(phi1)*exp(phi2); sigma[2,2] = exp(phi2)^2;
}

model{
  // 全てのパラメータの事前分布を独立な正規分布で仮定する
  alpha_0 ~ multi_normal(rep_vector(0,2),diag_matrix(rep_vector(1,2))); // ~N_2((0,0),(1,0,0,1) )
  for(i in 1:2*P){
    alpha_1[1,i] ~ normal(0,1); // N(0,1)
    alpha_1[2,i] ~ normal(0,1); // N(0,1)
  }
  phi1 ~ normal(0,10^5); phi2 ~ normal(0,10^5); phi3 ~ normal(0,10^5);
  for(n in 1+P:N){
    vector[P] pre_theta; // P期前までのtheta ベクトルを用意する.
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    theta[n] ~ circular_reg_lpdf(P,pre_theta,alpha_0,alpha_1,sigma);
  }
}

generated quantities{
  vector[N-P] log_likelihood;
  for(n in 1+P:N){
    vector[P] pre_theta; // P期前までのtheta ベクトルを用意する.
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    log_likelihood[n-P] = circular_reg_lpdf(theta[n]| P, pre_theta, alpha_0, alpha_1, sigma);
  } 
}
