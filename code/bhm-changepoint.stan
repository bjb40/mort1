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
  int<lower=0> YRS; //n year-specific effects
  int yrctr; //integer to center year for index 1: end-year
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
  vector<lower=0>[TDS] zi; // scale for correlation matrix
  vector<lower=0>[TDS] delta; //scale for year-specific errors
  cholesky_factor_corr[2] L_Omega; //faster for programming; correlation matrix
  matrix[TDS,IDS] omega_i; //container for random normal draw to distribute cross-cell error
  matrix[TDS,YRS] omega_t; //container for random normal draw to distribut cross-year error
}

transformed parameters {
    matrix[TDS,IDS] mu_i; // age-specific conditonal effects
    matrix[TDS,YRS] mu_t; // year-specific conditional effects
    vector[N] yhat;
    mu_i <- diag_matrix(zi)*L_Omega*omega_i;
    mu_t[1] <- delta[1]*omega_t[1];
    mu_t[2] <- delta[2]*omega_t[2];

for(n in 1:N){
    yhat[n] <- mu_i[td[n],id[n]] + mu_t[td[n],t[n]+yrctr] + t[n]*beta[td[n]] + z[n]*gamma[td[n]]';
  }
  
}

model{

  to_vector(omega_i) ~ normal(0,1);
  to_vector(omega_t) ~ normal(0,1);

  y ~ normal(yhat,sig);

  //prior
  beta ~ normal(0,5);
  to_vector(gamma) ~ normal(0,5);
  sig ~ normal(0,5);
  zi ~ cauchy(0,5);
  delta ~ cauchy(0,5);
  L_Omega ~ lkj_corr_cholesky(1); //1 is equiv to uniform prior; >1 diagonal <1 high

}

// see DIC in stan's google mailing list for discussion, BDA3, pp.172-179

generated quantities {
  //for WAIC
  vector[N] loglik; // log pointwise predictive density
  //for DIC
  real dev;

  dev <- 0;  
  for(n in 1:N){
    loglik[n] <- (normal_log(y[n],yhat[n],sig));
    dev <- dev-(2*normal_log(y[n],yhat[n],sig));

  }
}
