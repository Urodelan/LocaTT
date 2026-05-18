//// Multivariate logistic regression

functions{
  
  //// Define custom functions.
  
  // Define the log probability density function of the Gaussian copula.
  real copula_cholesky_lpdf(vector u, matrix L){
    vector[rows(u)] z;
    z = inv_Phi(u);
    return multi_normal_cholesky_lpdf(z|rep_vector(0,rows(z)),L) - std_normal_lpdf(z);
  }
  
  // Define the log probability density function of the multivariate logistic.
  real multi_logistic_cholesky_lpdf(vector y, vector mu, vector sigma, matrix L){
    real density;
    vector[rows(y)] u;
    density = 0;
    for(i in 1:rows(y)){
      density += logistic_lpdf(y[i]|mu[i],sigma[i]);
      u[i] = logistic_cdf(y[i]|mu[i],sigma[i]);
    }
    density += copula_cholesky_lpdf(u|L);
    return density;
  }
  
  // Define sum function for 2D integer array.
  int sum2d(array[,] int a){
    int s = 0;
    for(i in 1:size(a)){
      s += sum(a[i]);
    }
    return s;
  }
  
}

data{
  
  //// Define variables.
  int<lower=1> N; // Number of observations.
  int<lower=1> D; // Number of dimensions.
  int<lower=1> K; // Number of predictors.
  int<lower=0,upper=1> Y[N,D]; // Response matrix.
  matrix[N,K] X; // Predictor matrix.
  
  //// Define priors.
  real<lower=0> lkj; // LKJ correlation prior.
  real B_mu; // Regression coefficient mean prior.
  real<lower=0> B_sd; // Regression coefficient standard deviation prior.
  
}

transformed data {
  
  // Specify variables.
  int<lower=0> N_pos;
  array[sum2d(Y)] int<lower=1,upper=N> n_pos;
  array[size(n_pos)] int<lower=1,upper=D> d_pos;
  int<lower=0> N_neg;
  array[(N*D)-size(n_pos)] int<lower=1,upper=N> n_neg;
  array[size(n_neg)] int<lower=1,upper=D> d_neg;
  
  // Define variables.
  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i;
    int j;
    i = 1;
    j = 1;
    for(n in 1:N){
      for(d in 1:D){
        if(Y[n,d]==1){
          n_pos[i] = n;
          d_pos[i] = d;
          i += 1;
        }else{
          n_neg[j] = n;
          d_neg[j] = d;
          j += 1;
        }
      }
    }
  }
  
}

parameters{
  
  //// Specify parameters.
  matrix[D,K] B; // Regression coefficient matrix.
  vector<lower=0>[N_pos] z_pos; // Positive latent variable.
  vector<upper=0>[N_neg] z_neg; // Negative latent variable.
  cholesky_factor_corr[D] L; // Cholesky factor of correlation matrix.
  
}

transformed parameters{
  
  //// Specify transformed parameters.
  array[N] vector[D] z; // Latent variable.
  
  //// Define transformed parameters.
  
  // Latent variables.
  for(n in 1:N_pos){
    z[n_pos[n],d_pos[n]] = z_pos[n];
  }
  for(n in 1:N_neg){
    z[n_neg[n],d_neg[n]] = z_neg[n];
  }
  
}

model{
  
  //// Provide priors.
  
  // Regression coefficients.
  for(d in 1:D){
    for(k in 1:K){
      target += normal_lpdf(B[d,k]|B_mu,B_sd);
    }
  }
  
  // Cholesky factor of correlation matrix.
  target += lkj_corr_cholesky_lpdf(L|lkj);
  
  //// State likelihood.
  
  // Loop through each observation.
  for(i in 1:N){
    
    // Define local variable.
    vector[D] eta;
    
    // Loop through each dimension.
    for(d in 1:D){
      
      // Linear model for eta.
      eta[d] = B[d,1:K] * transpose(X[i,1:K]);
      
    }
    
    // Multivariate logistic likelihood.
    target += multi_logistic_cholesky_lpdf(z[i]|eta,rep_vector(1,D),L);
    
  }
  
}

generated quantities{
  
  //// Specify generated parameters.
  corr_matrix[D] R; // Residual correlation matrix.
  
  //// Define generated parameters.
  
  // Residual correlation matrix.
  R = multiply_lower_tri_self_transpose(L);
  
}
