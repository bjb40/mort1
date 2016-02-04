#@@@@@@@@@@@@@@@@@@@
#modified changepoint model/simplifed interaction
#this is really a multigroup hierarchical, b/c we know the changepoint
#but iti s like the Carlin 1992 model on p. 401, with k known
#dev Stan 2.7.0
#Bryce Bartlett
#@@@@@@@@@@@@@@@@@@@

data { 
  int<lower=0> N; //n observations
  int<lower=0> IDS; //n cells
  int<lower=0> TDS; //n groups (2) in this case
  int id[N]; // unique age groups (for random effect)
  real y[N]; // outcomes
  int t[N]; //time variable (years)
  int td[N]; //index number for gorup of time
  int<lower=0> P; //dimensions of predictors
  matrix[N,P] z; // all time invariant
  } 
  
parameters{
  #individual level
  vector[2] beta; // grand mean coefficients for intercept and slope
  matrix[TDS,P] gamma; //
  real<lower=0> sig;//l1 error; BDA3 388 - uniform gelman 2006; stan manual 66
  vector<lower=0>[TDS]zi; // scale for correlation matrix
  cholesky_factor_corr[2] L_Omega; //faster for programming; correlation matrix
  matrix[TDS,IDS] omega_i; //container for random normal draw to distribute cross-cell error
}

transformed parameters {
    matrix[TDS,IDS] mu_i; // age-specific conditonal effects
    vector[N] yhat;
    mu_i <- diag_matrix(zi)*L_Omega*omega_i;
    

  for(n in 1:N){
    yhat[n] <- mu_i[td[n],id[n]] + t[n]*beta[td[n]] + z[n]*gamma[td[n]]';
  }

}

model{

  to_vector(omega_i) ~ normal(0,1);

  y ~ normal(yhat,sig);
  
  //prior
  beta ~ normal(0,5);
  to_vector(gamma[1]) ~ normal(0,5);
  to_vector(gamma[2]) ~ normal(gamma[2],5);
  
  sig ~ normal(0,5);

  zi ~ cauchy(0,5);

}

// see DIC in stan's google mailing list for discussion, BDA3, pp.172-179

generated quantities {
  //for WAIC
  vector[N] loglik; // log pointwise predictive density
  //for DIC
  real dev;
  //FOR PPD
  vector[N] ppd;
  matrix[2,2] sigma; //covariance matrix
  
  sigma <- quad_form_diag(L_Omega*L_Omega',zi);

  dev <- 0;  
  for(n in 1:N){
    loglik[n] <- (normal_log(y[n],yhat[n],sig));
    dev <- dev-(2*normal_log(y[n],yhat[n],sig));
    ppd[n] <- normal_rng(yhat[n],sig);
  }
}
