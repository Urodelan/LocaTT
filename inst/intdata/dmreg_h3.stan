//// Dirichlet-multinomial regression
// Hierarchical effects: 3

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
  int<lower=1> L1; // Number of hierarchical effects (#1).
  int<lower=1> h1[N]; // Vector of hierarchical identifiers (#1).
  int<lower=1> L2; // Number of hierarchical effects (#2).
  int<lower=1> h2[N]; // Vector of hierarchical identifiers (#2).
  int<lower=1> L3; // Number of hierarchical effects (#3).
  int<lower=1> h3[N]; // Vector of hierarchical identifiers (#3).
  
  //// Define priors.
  real B_mu; // Fixed effects prior mean.
  real<lower=0> B_sd; // Fixed effects prior standard deviation.
  real theta_mu; // Precision parameter prior mean.
  real<lower=0> theta_sd; // Precision parameter prior standard deviation.
  real<lower=0> sigma2_alpha; // Hierarchical variance alpha prior.
  real<lower=0> sigma2_beta; // Hierarchical variance beta prior.
  
}

parameters{
  
  //// Specify parameters.
  matrix[D-1,K] B_raw; // Raw fixed effects matrix.
  real theta; // Precision parameter.
  matrix[D-1,L1] H1_raw; // Raw hierarchical effects matrix (#1).
  vector<lower=0>[D-1] sigma2_h1; // Hierarchical variance vector (#1).
  matrix[D-1,L2] H2_raw; // Raw hierarchical effects matrix (#2).
  vector<lower=0>[D-1] sigma2_h2; // Hierarchical variance vector (#2).
  matrix[D-1,L3] H3_raw; // Raw hierarchical effects matrix (#3).
  vector<lower=0>[D-1] sigma2_h3; // Hierarchical variance vector (#3).
  
}

transformed parameters{
  
  //// Specify transformed parameters.
  matrix[D,K] B; // Fixed effects matrix.
  real exptheta; // Exponentiated precision parameter.
  matrix[D,L1] H1; // Hierarchical effects matrix (#1).
  vector<lower=0>[D-1] sigma1; // Hierarchical standard deviation vector (#1).
  matrix[D,L2] H2; // Hierarchical effects matrix (#2).
  vector<lower=0>[D-1] sigma2; // Hierarchical standard deviation vector (#2).
  matrix[D,L3] H3; // Hierarchical effects matrix (#3).
  vector<lower=0>[D-1] sigma3; // Hierarchical standard deviation vector (#3).
  
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
  
  // Hierarchical effects (#1).
  // Set last dimension's effects to zero.
  for(l in 1:L1){
    H1[D,l]=0;
  }
  // Set remaining effects to raw coefficients.
  for(d in 1:(D-1)){
    for(l in 1:L1){
      H1[d,l]=H1_raw[d,l];
    }
  }
  // Calculate standard deviations from variances.
  sigma1=sqrt(sigma2_h1);
  
  // Hierarchical effects (#2).
  // Set last dimension's effects to zero.
  for(l in 1:L2){
    H2[D,l]=0;
  }
  // Set remaining effects to raw coefficients.
  for(d in 1:(D-1)){
    for(l in 1:L2){
      H2[d,l]=H2_raw[d,l];
    }
  }
  // Calculate standard deviations from variances.
  sigma2=sqrt(sigma2_h2);
  
  // Hierarchical effects (#3).
  // Set last dimension's effects to zero.
  for(l in 1:L3){
    H3[D,l]=0;
  }
  // Set remaining effects to raw coefficients.
  for(d in 1:(D-1)){
    for(l in 1:L3){
      H3[d,l]=H3_raw[d,l];
    }
  }
  // Calculate standard deviations from variances.
  sigma3=sqrt(sigma2_h3);
  
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
  
  // Hierarchical effects (#1).
  for(d in 1:(D-1)){
    for(l in 1:L1){
      target += normal_lpdf(H1_raw[d,l]|0,sigma1[d]);
    }
  }
  
  // Hierarchical variances (#1).
  for(d in 1:(D-1)){
    target += inv_gamma_lpdf(sigma2_h1[d]|sigma2_alpha,sigma2_beta);
  }
  
  // Hierarchical effects (#2).
  for(d in 1:(D-1)){
    for(l in 1:L2){
      target += normal_lpdf(H2_raw[d,l]|0,sigma2[d]);
    }
  }
  
  // Hierarchical variances (#2).
  for(d in 1:(D-1)){
    target += inv_gamma_lpdf(sigma2_h2[d]|sigma2_alpha,sigma2_beta);
  }
  
  // Hierarchical effects (#3).
  for(d in 1:(D-1)){
    for(l in 1:L3){
      target += normal_lpdf(H3_raw[d,l]|0,sigma3[d]);
    }
  }
  
  // Hierarchical variances (#3).
  for(d in 1:(D-1)){
    target += inv_gamma_lpdf(sigma2_h3[d]|sigma2_alpha,sigma2_beta);
  }
  
  //// State likelihood.
  
  // Loop through each observation.
  for(i in 1:N){
    
    // Define local variable.
    vector[D] eta;
    
    // Loop through each dimension.
    for(d in 1:D){
      // Linear model for eta.
      eta[d]=B[d,1:K]*transpose(X[i,1:K])+H1[d,h1[i]]+H2[d,h2[i]]+H3[d,h3[i]];
    }
    
    // Dirichlet-multinomial likelihood.
    target += dirichlet_multinomial_lpmf(Y[i,1:D]|softmax(eta[1:D])*exptheta);
    
  }
  
}
