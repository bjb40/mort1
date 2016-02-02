#@@@@@@@@@@@@@@@@@@@
#basic hierarchical model (models 1 and 3)
#dev Stan 2.7.0
#Bryce Bartlett
#@@@@@@@@@@@@@@@@@@@

data { 
  int<lower=0> N; //n observations
  int<lower=0> IDS; //n cells
  int id[N]; // unique age groups (for random effect)
  real y[N]; // outcomes
  int t[N]; //time variable (years)
  int<lower=0> P; //dimensions of predictors
  matrix[N,P] z; // all time invariant
  } 

parameters{
  #individual level
  vector[2] beta; // grand mean coefficients for intercept and slope
  vector[P] gamma; //
  real<lower=0> sig; //l1 error; BDA3 388 - uniform gelman 2006; stan manual 66
  real<lower=0> zi; //(scale for intercept)
  vector[IDS] omega_i; //container for random normal draw to distribute cross-cell error
}

transformed parameters {
    vector[IDS] mu_i; // age-specific conditonal effects
    mu_i <- beta[1] + omega_i*zi;

}

model{
  real yhat[N];
  
  to_vector(omega_i) ~ normal(0,1);

  for(n in 1:N){
    yhat[n] <- mu_i[id[n]] + t[n]*beta[2] + z[n]*gamma;
  }
  
    y ~ normal(yhat,sig);
  
  //prior
  to_vector(beta) ~ normal(0,5);
  gamma ~ normal(0,5);
  sig ~ normal(0,5);

  zi ~ cauchy(0,10);
  
}

// see DIC in stan's google mailing list for discussion, BDA3, pp.172-179

#generated quantities {
  //for WAIC
  #real PoinLikelihoods;
  #vector[N] PointLikelihoods; // log pointwise predictive density
  #for(i in 1:N){
  #  PointLikelioods[i] <- exp(normal_log(yhat[i],mu,sig))
  #}
  //for DIC
  #real dev;
  #dev <- 0;
  #for (i in 1:N){
  #  dev <- dev-(2*normal_log(y[i],yhat[i],sig))
  #}
  
#}
