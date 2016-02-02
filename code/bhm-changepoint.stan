#@@@@@@@@@@@@@@@@@@@
#modified changepoint model/simplifed interaction
#dev Stan 2.7.0
#Bryce Bartlett
#@@@@@@@@@@@@@@@@@@@

data { 
  int<lower=0> N; //n observations
  int<lower=0> IDS; //n cells
  int<lower=0> TDS; //n years
  int id[N]; // unique age groups (for random effect)
  real y[N]; // outcomes
  int t[N]; //time variable (years)
  int td[N]; //index number of time variable (random effect)
  int<lower=0> P; //dimensions of predictors
  matrix[N,P] z; // all time invariant
  } 

transformed data{
  vector[IDS] i_int;
  i_int <- rep_vector(1.0,IDS);
}

parameters{
  #individual level
  vector[2] beta; // grand mean coefficients for random effects
  vector[P] gamma; //
  real<lower=0> sig; //l1 error; BDA3 388 - uniform gelman 2006; stan manual 66
  vector<lower=0>[2] zi; //(scale for intercept and slope
  vector[IDS] omega_i; //container for random normal draw to distribute age error
  vector[TDS] omega_t; //container for random normal draw to distirbute time error
}

transformed parameters {
    vector[IDS] mu_i; // age-specific conditonal effects
    vector[TDS] mu_t; // time-specific conditional effects
    mu_i <- i_int*beta[1] + omega_i*zi[1];
    mu_t <- omega_t*zi[2];
  
}

model{
  real yhat[N];
  
  to_vector(omega_i) ~ normal(0,1);
  to_vector(omega_t) ~ normal(0,1);
  
  for(n in 1:N){
    yhat[n] <- mu_i[id[n]] + t[n]*beta[2] + z[n]*gamma + mu_t[td[n]];
  }
  
    y ~ normal(yhat,sig);
  
  //prior
  to_vector(beta) ~ normal(0,5);
  gamma ~ normal(0,5);
  sig ~ normal(0,5);

  zi ~ normal(0,1.5);
  
}
