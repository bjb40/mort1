#@@@@@@@@@@@@@@@@@@@
#basic hierarchical model (models 1 and 3)
#dev Stan 2.7.0
#Bryce Bartlett
#@@@@@@@@@@@@@@@@@@@

data { 
  int<lower=0> N; //n observations
  int<lower=0> IDS; //n cells
  int<lower=0> YRS; //n years
  int yrctr; //number to add to t to create an index of integers from 1 to end
  int id[N]; // unique age groups (for random effect)
  real y[N]; // outcomes
  int t[N]; //time variable (years)
  int<lower=0> P; //dimensions of predictors
  matrix[N,P] z; // all time invariant
  } 

parameters{
  #individual level
  real beta; // grand mean coefficient for slope
  vector[P] gamma; //
  real<lower=0> sig; //l1 error; BDA3 388 - uniform gelman 2006; stan manual 66
  real<lower=0> zi; //(scale for intercept)
  real<lower=0> delta; // variance for years
  vector[IDS] omega_i; //container for random normal draw to distribute cross-cell error
  vector[YRS] omega_t; //container for random normal draw across year
}

transformed parameters {
    vector[IDS] mu_i; // age-specific conditonal effects
    vector[YRS] mu_t; // year-specific effects
    vector[N] yhat;
    mu_i <- omega_i*zi;
    mu_t <- omega_t*delta;

  for(n in 1:N){
    yhat[n] <- mu_i[id[n]] + mu_t[t[n] +yrctr] + t[n]*beta + z[n]*gamma;
  }

}

model{

  to_vector(omega_i) ~ normal(0,1);
  to_vector(omega_t) ~ normal(0,1); // can overdisperse as necessary, but may need reparameter

  
    y ~ normal(yhat,sig);
  
  //prior
  beta ~ normal(0,5);
  gamma ~ normal(0,5);
  sig ~ normal(0,5);

  zi ~ cauchy(0,5);
  delta ~ cauchy(0,5);
  
}

// see DIC in stan's google mailing list for discussion, BDA3, pp.172-179

generated quantities {
  //for WAIC
  vector[N] loglik; // log pointwise predictive density
  //for DIC
  real dev;
  //FOR PPD
  vector[N] ppd;

  dev <- 0;  
  for(i in 1:N){
    loglik[i] <- (normal_log(y[i],yhat[i],sig));
    dev <- dev-(2*normal_log(y[i],yhat[i],sig));
    ppd[i] <- normal_rng(yhat[i],sig);
  }

}
