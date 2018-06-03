functions{
  real circular_reg_lpdf(real theta, int P, vector mu_hat, matrix sigma){
    vector[2] u; real A; real B; real C; real D; real p;
    u[1] = cos(theta); u[2] = sin(theta);
    A = quad_form(inverse_spd(sigma), u); B = u' * inverse_spd(sigma) * mu_hat;
    C = (-0.5) * quad_form(inverse_spd(sigma), mu_hat); D = B/sqrt(A);
    p = -log(A) - 0.5*log(determinant(sigma)) + C
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
  vector[2] alpha_0;
  matrix[2,2*P] alpha_1; // P個の係数行列
  real phi1;
  real phi2;
  real phi3;
}

transformed parameters{
  unit_vector[2] mu_hat[N-P];
  cov_matrix[2] sigma;
  vector[2] mu;
  // Define mu_hat 
  for(n in 1+P:N){
    vector[2*P] tmp;
    for(k in 1:P){
      tmp[2*k - 1] = cos(theta[n-k]); 
      tmp[2*k] = sin(theta[n-k]);
    }
    mu = alpha_0 + ( alpha_1 * tmp);
    mu_hat[n-P] = mu / sqrt( mu' * mu);
    // print("tmp=",tmp, ",mu=",mu, ",mu_hat=",mu_hat[n-P]); //debug
  }
  // Define Sigma
  sigma[1,1] = exp(phi1)^2; sigma[1,2] = tanh(phi3)*exp(phi1)*exp(phi2);
  sigma[2,1] = tanh(phi3)*exp(phi1)*exp(phi2); sigma[2,2] = exp(phi2)^2;
}

model{
  // 全てのパラメータの事前分布を独立な正規分布で仮定する
  alpha_0 ~ multi_normal(rep_vector(0,2),diag_matrix(rep_vector(10^5,2))); // ~N_2((0,0),(10^5,0,0,10^5) )
  for(i in 1:2*P){
    alpha_1[1,i] ~ normal(0,10^5); // N(0,10^5)
    alpha_1[2,i] ~ normal(0,10^5); // N(0,10^5)
  }
  phi1 ~ normal(0,10^5); phi2 ~ normal(0,10^5); phi3 ~ normal(0,10^5);
  // 尤度を最大化する
  for(n in 1+P:N){
    theta[n] ~ circular_reg_lpdf(P,mu_hat[n-P], sigma);
  }
}

// generated quantities{
//   vector[N-P] log_likelihood;
//   for(n in 1+P:N){
//     vector[P] pre_theta; // P期前までのtheta ベクトルを用意する.
//     for(k in 1:P){
//       pre_theta[k] = theta[n-k];
//     }
//     log_likelihood[n-P] = circular_reg_lpdf(theta[n]| P, pre_theta, alpha_0, alpha_1, sigma);
//   } 
// }

