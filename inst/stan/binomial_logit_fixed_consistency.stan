// Fixed effects consistency model for binomial data with logit link
// Network meta-analysis with optional meta-regression

data {
  int<lower=0> ns_a;                      // number of arm-based studies
  int<lower=0> nt;                        // number of treatments
  array[ns_a] int<lower=2> na_a;          // number of arms per study
  array[ns_a, max(na_a)] int<lower=0> r;  // number of events
  array[ns_a, max(na_a)] int<lower=0> n;  // number of participants
  array[ns_a, max(na_a)] int<lower=1, upper=nt> t_a;  // treatment indices
  
  // Priors
  real prior_mu_mean;
  real<lower=0> prior_mu_sd;
  real prior_d_mean;
  real<lower=0> prior_d_sd;
  
  // Meta-regression (optional)
  int<lower=0, upper=1> has_metareg;      // indicator for meta-regression
  array[ns_a, max(na_a)] real x_a;        // covariate values (if has_metareg=1)
  
  int<lower=0, upper=3> metareg_type;     // 0=none, 1=UNRELATED, 2=EXCHANGEABLE, 3=EQUAL
  real prior_beta_mean;
  real<lower=0> prior_beta_sd;
}

transformed data {
  int max_arms = max(na_a);
}

parameters {
  vector[ns_a] mu;                        // baseline effects
  vector[nt-1] d_raw;                     // treatment effects (relative to reference)
  
  // Meta-regression parameters
  vector[nt-1] beta_raw;                  // regression coefficients (if UNRELATED)
  real b;                                 // mean coefficient (if EXCHANGEABLE)
  real<lower=0> gamma;                    // sd of coefficients (if EXCHANGEABLE)
  real B;                                 // common coefficient (if EQUAL)
}

transformed parameters {
  vector[nt] d;
  vector[nt] beta;
  
  d[1] = 0;  // reference treatment
  d[2:nt] = d_raw;
  
  // Meta-regression coefficients
  beta[1] = 0;
  if (has_metareg == 0) {
    beta[2:nt] = rep_vector(0, nt-1);
  } else if (metareg_type == 1) {  // UNRELATED
    beta[2:nt] = beta_raw;
  } else if (metareg_type == 2) {  // EXCHANGEABLE
    beta[2:nt] = b + gamma * beta_raw;
  } else if (metareg_type == 3) {  // EQUAL
    beta[2:nt] = rep_vector(B, nt-1);
  }
}

model {
  // Priors
  mu ~ normal(prior_mu_mean, prior_mu_sd);
  d_raw ~ normal(prior_d_mean, prior_d_sd);
  
  // Meta-regression priors
  if (has_metareg == 1) {
    if (metareg_type == 1) {  // UNRELATED
      beta_raw ~ student_t(1, prior_beta_mean, prior_beta_sd);
    } else if (metareg_type == 2) {  // EXCHANGEABLE
      beta_raw ~ std_normal();
      b ~ student_t(1, prior_beta_mean, prior_beta_sd);
      gamma ~ uniform(0, prior_beta_sd);
    } else if (metareg_type == 3) {  // EQUAL
      B ~ student_t(1, prior_beta_mean, prior_beta_sd);
    }
  }
  
  // Likelihood
  for (i in 1:ns_a) {
    for (k in 1:na_a[i]) {
      real theta;
      theta = mu[i] + d[t_a[i,k]] - d[t_a[i,1]];
      if (has_metareg == 1) {
        theta += (beta[t_a[i,k]] - beta[t_a[i,1]]) * x_a[i,k];
      }
      r[i,k] ~ binomial_logit(n[i,k], theta);
    }
  }
}

generated quantities {
  array[ns_a, max_arms] real rhat;
  array[ns_a, max_arms] real dev_a;
  real totresdev = 0;
  
  // Calculate deviance contributions
  for (i in 1:ns_a) {
    for (k in 1:na_a[i]) {
      real theta;
      real p;
      theta = mu[i] + d[t_a[i,k]] - d[t_a[i,1]];
      if (has_metareg == 1) {
        theta += (beta[t_a[i,k]] - beta[t_a[i,1]]) * x_a[i,k];
      }
      p = inv_logit(theta);
      rhat[i,k] = p * n[i,k];
      
      if (r[i,k] > 0 && r[i,k] < n[i,k]) {
        dev_a[i,k] = 2 * (r[i,k] * (log(r[i,k]) - log(rhat[i,k])) + 
                          (n[i,k] - r[i,k]) * (log(n[i,k] - r[i,k]) - log(n[i,k] - rhat[i,k])));
      } else if (r[i,k] == 0) {
        dev_a[i,k] = 2 * (n[i,k] - r[i,k]) * (log(n[i,k] - r[i,k]) - log(n[i,k] - rhat[i,k]));
      } else {  // r[i,k] == n[i,k]
        dev_a[i,k] = 2 * r[i,k] * (log(r[i,k]) - log(rhat[i,k]));
      }
      totresdev += dev_a[i,k];
    }
  }
}

