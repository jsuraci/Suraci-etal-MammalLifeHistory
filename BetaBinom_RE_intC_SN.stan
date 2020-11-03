//-----------------------------------------------------------------------
//          ** BETA-BINOMIAL OCCUPANCY MODEL SS BS GROUP1 **
// SINGLE SPECIES - RANDOM EFFECTS MODEL -  WITH CAM RE
//-----------------------------------------------------------------------
// DESCRIPTION: Single-species occupancy model with project 
// (i.e., group) level coefficients for the effects of human 
// presence and human footprint on Occupancy Probability and detection. 
// Project-level coefficients are RANDOM (i.e., unmodeled) 
// Camera-level random intercepts INCLUDED and are nested in project
// Posterior predictive distribution (y_new) and stats for bayesian 
// p-values are calculated outside of the model so excluded here.
// Formulated as **BETA-BINOMIAL** occupancy model with rho 
// overdispersion parameter

data {
  
  int<lower = 0> N; // Number of observations (i.e., sampling periods at camera sites)
  int<lower = 1> J; // Number of unique projects
  int<lower = 1> C; // Number of camera sites
  int<lower = 1> Kp; // Number of observation-level covariates WITH PROJECT COEFS
  int<lower = 1> Kfa; // Number of FIXED observation-level covariates on occupancy
  int<lower = 1> Kfb; // Number of FIXED observation-level covariates on detection
  int<lower=1,upper=J> projO[N]; // Codes which project each observation in N belongs to
  int<lower = 1,upper = J> projC[C]; // Codes which project each camera site in C belongs to
  int<lower = 1, upper = C> camO[N]; // Codes which camera each observation in N belongs to

  int<lower = 0> y[N]; // Detection data (all independent visits)
  real<lower=0, upper=1> Q[N]; // Naive occupancy (at least one detection during S weeks)
  int<lower = 1> S[N]; // Number of weeks (binomial trials)

  // Observation covariates
  row_vector[Kp] Xp[N]; // Observation covariate data with PROJECT-LEVEL COEFS (on occupancy and detection)
  row_vector[Kfa] Xfa[N]; // FIXED observation covariate data on occupancy
  row_vector[Kfb] Xfb[N]; // FIXED observation covariate data on detection

}

parameters {
  
  // Define project-level coefficients as raw values 
  // and their group-level means variances
  matrix[J, Kp] A_raw;
  matrix[J, Kp] B_raw;
  real A_mu[Kp];
  real B_mu[Kp];
  real <lower = 0> A_sigma[Kp];
  real <lower = 0> B_sigma[Kp];
  
  // Define fixed coefficients on observation data
  vector[Kfa] a;
  vector[Kfb] b;
  
   // Define camera site random intercept plus variance of dist
  vector[C] a0_raw; // camera-level random intercept for occupancy (raw value for SN param)
  vector[C] b0_raw; // camera-level random intercept for detection (raw value for SN param)
  real <lower = 0> sig_a0; // Variance for distribution for cam level intercepts
  real <lower = 0> sig_b0; // Variance for distribution for cam level intercepts
  
  // Define overdispersion paremeter rho
  real<lower = 0, upper = 1> rho;
  
}

transformed parameters {
  
  // Define occupancy probability (psi) and poisson parameters for visitation
  vector<lower = 0, upper = 1>[N] psi;
  vector<lower = 0, upper = 1>[N] mu;
  vector<lower = 0>[N] alpha;
  vector<lower = 0>[N] beta;
  
  // Define true values for parameters using SN parameterization
  matrix[J, Kp] A;
  matrix[J, Kp] B;
  vector[C] a0; 
  vector[C] b0;
  
  // Create project-level coefficients on observation data (true values) from raw values
  for(k in 1:Kp){
    A[,k] = A_mu[k] + (A_sigma[k] * A_raw[,k]);
    B[,k] = B_mu[k] + (B_sigma[k] * B_raw[,k]);
  }
  
  // Create camera site random intercepts (true values) from raw values
  // Mean of cam level intercept is project level intercept
  a0 = A[projC,1] + (sig_a0 * a0_raw);
  b0 = B[projC,1] + (sig_b0 * b0_raw);
  
  for(i in 1:N){
    // Model psi and p as a function of observation-level covariates
    mu[i] = inv_logit(b0[camO[i]] + dot_product(B[projO[i], 2:Kp], Xp[i, 2:Kp]) + dot_product(b, Xfb[i]));
    psi[i] = inv_logit(a0[camO[i]] + dot_product(A[projO[i], 2:Kp], Xp[i, 2:Kp]) + dot_product(a, Xfa[i]));
    
    // Calculate parameters of beta-binomial dist from mu and rho
    alpha[i] = mu[i] * (1 - rho) / rho ;
    beta[i] = (1 - mu[i]) * (1 - rho) / rho; 
  }

}

model {
  
  // Priors
  a ~ normal(0, 1); // Vectorized over Kfa
  b ~ normal(0, 1); // Vectorized over Kfb
  sig_a0 ~ normal(0, 1);
  sig_b0 ~ normal(0, 1);
  a0_raw ~ std_normal(); 
  b0_raw ~ std_normal();
  
  A_mu ~ normal(0, 1); // Get priors and values for A and B - vectorized over K
  B_mu ~ normal(0, 1);
  A_sigma ~ normal(0, 1);
  B_sigma ~ normal(0, 1);
  for (j in 1:J){
      A_raw[j] ~ std_normal(); // Standard normal priors for A_raw and B_raw
      B_raw[j] ~ std_normal();
  }
  
  rho ~ uniform(0,1); // Prior for overdispersion parameter

  
  // Likelihood
  for(i in 1:N){
    if(Q[i] == 1){
      target += bernoulli_lpmf(1 | psi[i])
                + beta_binomial_lpmf(y[i] | S[i], alpha[i], beta[i]);
    }
    if(Q[i] == 0){
      target += log_sum_exp(bernoulli_lpmf(1 | psi[i]) + beta_binomial_lpmf(0 | S[i], alpha[i], beta[i]), bernoulli_lpmf(0 | psi[i]));
    }
  }
}

generated quantities{
  // Declare new variables
  vector[N] log_lik; // Model log likelihood
    
  // Get pointwise model log likelihood
  for(i in 1:N){
    if(Q[i] == 1){
      log_lik[i] = bernoulli_lpmf(1 | psi[i])
                + beta_binomial_lpmf(y[i] | S[i], alpha[i], beta[i]);
    }
    if(Q[i] == 0){
      log_lik[i] = log_sum_exp(bernoulli_lpmf(1 | psi[i]) + beta_binomial_lpmf(0 | S[i], alpha[i], beta[i]), bernoulli_lpmf(0 | psi[i]));
    }
  }
}
