//// Dirichlet-multinomial regression
// Hierarchical effects: 0

functions{
  
  //// Define custom functions.
  
  // Define Dirichlet-multinomial log probability mass function.
  real dirichlet_multinomial_lpmf(int[] y, vector alpha){
    return lgamma(sum(alpha))+lgamma(sum(y)+1)-lgamma(sum(y)+sum(alpha))+sum(lgamma(to_vector(y)+alpha))-sum(lgamma(alpha))-sum(lgamma(to_vector(y)+1));
  }
  
}

data{
  
  //// Define variables.
  int<lower=1> N; // Number of observations.
  int<lower=1> D; // Number of dimensions.
  int<lower=1> K; // Number of predictors.
  int Y[N,D]; // Response matrix.
  matrix[N,K] X; // Predictor matrix.
  
  //// Define priors.
  real B_mu; // Fixed effects prior mean.
  real<lower=0> B_sd; // Fixed effects prior standard deviation.
  real theta_mu; // Precision parameter prior mean.
  real<lower=0> theta_sd; // Precision parameter prior standard deviation.
  
}

parameters{
  
  //// Specify parameters.
  matrix[D-1,K] B_raw; // Raw fixed effects matrix.
  real theta; // Precision parameter.
  
}

transformed parameters{
  
  //// Specify transformed parameters.
  matrix[D,K] B; // Fixed effects matrix.
  real exptheta; // Exponentiated precision parameter.
  
  //// Define transformed parameters.
  
  // Fixed effects.
  for(k in 1:K){
    // Set effects to raw coefficients.
    for(d in 1:(D-1)){
      B[d,k]=B_raw[d,k];
    }
    // Set last dimension's effects to zero.
    B[D,k]=0;
  }
  
  // Exponentiate precision parameter.
  exptheta=exp(theta);
  
}

model{
  
  //// Provide priors.
  
  // Fixed effects.
  for(k in 1:K){
    for(d in 1:(D-1)){
      target += normal_lpdf(B_raw[d,k]|B_mu,B_sd);
    }
  }
  
  // Precision parameter.
  target += normal_lpdf(theta|theta_mu,theta_sd);
  
  //// State likelihood.
  
  // Loop through each observation.
  for(i in 1:N){
    
    // Define local variable.
    vector[D] eta;
    
    // Loop through each dimension.
    for(d in 1:D){
      // Linear model for eta.
      eta[d]=B[d,1:K]*transpose(X[i,1:K]);
    }
    
    // Dirichlet-multinomial likelihood.
    target += dirichlet_multinomial_lpmf(Y[i,1:D]|softmax(eta[1:D])*exptheta);
    
  }
  
}
