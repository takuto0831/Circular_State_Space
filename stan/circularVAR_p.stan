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
    A = quad_form(inverse_spd(sigma), u); B = u' * inverse(sigma) * mu;
    C = (-0.5) * quad_form(inverse_spd(sigma), mu); D = B/sqrt(A);
    p = -log(A) - 0.5*log(determinant(sigma)) + C
    + log(1+(D * normal_cdf(D,0,1)/exp(normal_lpdf(D|0,1))));    
    return p;
  }
}
  
data{
  int N; // sample size
  int P; // VAR(P) 
  real<lower=0,upper=2*pi()> theta[N]; // data
}

parameters{
  unit_vector[2] alpha_0;
  matrix[2,2*P] alpha_1; // P個の係数行列
  cov_matrix[2] sigma;
}

model{
  alpha_0 ~ multi_normal(rep_vector(0,2),diag_matrix(rep_vector(10^5,2)));
  sigma ~ inv_wishart(2,diag_matrix(rep_vector(1,2)));
  for(n in 1+P:N){
    vector[P] pre_theta; // P期前までのtheta ベクトルを用意する.
    for(k in 1:P){
      pre_theta[k] = theta[n-k];
    }
    theta[n] ~ circular_reg_lpdf(P,pre_theta,alpha_0,alpha_1,sigma);
  }
}

