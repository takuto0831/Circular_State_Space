functions{
  real circular_reg_lpdf(real theta, vector mu, matrix sigma){
    vector[2] u;
    real A;
    real B;
    real C;
    real D;
    real p;
    u[1] = cos(theta); u[2] = sin(theta);
    A = u' * inverse(sigma) * u; B = u' * inverse(sigma) * mu;
    C = (-0.5) * (mu' * inverse(sigma) * mu); D = B/sqrt(A);
    p = -log(A)-log(sqrt(determinant(sigma))) + C
    + log(1+(D * normal_cdf(D,0,1)/exp(normal_lpdf(D|0,1))));    
    return p;
  }
}
  
  data{
  int N; // sample size
  int P; // lag
  real<lower=0,upper=2*pi()> theta[N]; //data
  }
  
  parameters{
  vector[2] mu;
  real<lower=0.0001> tau;
  real rho;
  }
