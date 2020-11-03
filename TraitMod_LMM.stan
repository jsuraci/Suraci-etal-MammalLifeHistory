//-----------------------------------------------------------------------
//                      ** TRAIT RESPONSE MODEL **
//-----------------------------------------------------------------------


data {
  
  int<lower = 0> N; // Number of observations (i.e., species x projects)
  int<lower = 1> J; // Number projects
  int<lower = 1> S; // Number of species
  int<lower = 1> Fa; // Number of families
  int Kf; // Number of covariates with fixed coefficients
  int<lower=1,upper=S> spO[N]; // Codes which species each observation in N belongs to
  int<lower=1,upper=Fa> famO[N]; // Codes which family each observation in N belongs to
  int famSp[S]; // Codes which family each species belongs to
  real y[N]; // Response data
  real<lower = 0> sigma_o[N]; // Estimation variance for each species
  row_vector[Kf] Xf[N]; // Covariates with fixed coefs

}

parameters {
  
  // Define family-level random intercept as raw values. Mu comes from order 
  real Bf_mu;
  vector[Fa] Bf_raw;
  real <lower = 0> Bf_sigma;
  
  // Define species-level random intercept as raw values.  Mu comes from family
  vector[S] Bs_raw; // species-level random intercept for occupancy (raw value for SN param)
  real <lower = 0> Bs_sigma; // Variance for distribution for species level intercepts
  
  // Define fixed coefs
  vector[Kf] b;
  
  // Define the process and observation error terms
  real<lower = 0> sigma_p;
  
  // Declare "raw" value for z ("true" outcomes of process)
  real z_raw[N];
  
}

transformed parameters {
  // Declare new variables that are functions of raw values
  real z[N];
  vector[Fa] Bf;
  vector[S] Bs;
  real x_beta[N];
  
  // Calculate family- and species-level intercepts from raw values
  Bf = Bf_mu + (Bf_sigma * Bf_raw);
  Bs = Bf[famSp] + (Bs_sigma * Bs_raw);
  
  // Calculate "true" outcomes z from raw values
  for (i in 1:N){
    x_beta[i] = Bs[spO[i]] + dot_product(b, Xf[i]);
    z[i] = x_beta[i] + (sigma_p * z_raw[i]);
  }

}

model {

  // Priors
  z_raw ~ std_normal();
  Bf_mu ~ normal(0, 1000); // Get priors and values for A and B - vectorized over Kr
  Bf_sigma ~ uniform(0,100);
  Bf_raw ~ std_normal();
  Bs_sigma ~ uniform(0,100);
  Bs_raw ~ std_normal();
  b ~ normal(0, 1000); // Fixed interaction coef
  sigma_p ~ uniform(0,100); // Prior for process error
  
  // Likelihood (with efficiency trick from stan user guide, section 1.13)
  y ~ normal(z, sigma_o);
  
}

generated quantities{
  // Declare new variables
  vector[N] y_new;
  real cvar_new;
  vector[N] log_lik; // Model log likelihood

  // Simulate new data
  {
    real z_new;
    for(i in 1:N){
      z_new = normal_rng(x_beta[i], sigma_p);
      y_new[i] = normal_rng(z_new, sigma_o[i]);
      log_lik[i]  = normal_lpdf(y[i] | z[i], sigma_o[i]);
    }    
  }

  // Calcuate new data coefficient of variation
  cvar_new = sd(y_new)/mean(y_new);

}

