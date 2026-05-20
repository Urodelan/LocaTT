//// Univariate logistic regressions

data{
  
  //// Define variables.
  int<lower=1> N; // Number of observations.
  int<lower=1> D; // Number of dimensions.
  int<lower=1> K; // Number of predictors.
  int<lower=0,upper=1> Y[N,D]; // Response matrix.
  matrix[N,K] X; // Predictor matrix.
  
  //// Define priors.
  real B_mu; // Regression coefficient mean prior.
  real<lower=0> B_sd; // Regression coefficient standard deviation prior.
  
}

parameters{
  
  //// Specify parameters.
  matrix[D,K] B; // Regression coefficient matrix.
  
}

model{
  
  //// Provide priors.
  
  // Regression coefficients.
  for(d in 1:D){
    for(k in 1:K){
      target += normal_lpdf(B[d,k]|B_mu,B_sd);
    }
  }
  
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
    
    // Bernoulli likelihood (with logit link).
    target += bernoulli_logit_lpmf(Y[i]|eta);
    
  }
  
}
