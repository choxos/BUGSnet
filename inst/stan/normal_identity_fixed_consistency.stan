// Fixed effects consistency model for normal data with identity link
// For continuous outcomes (mean differences)

data {
  int<lower=0> ns_a;
  int<lower=0> nt;
  array[ns_a] int<lower=2> na_a;
  array[ns_a, max(na_a)] real y;  // outcome means
  array[ns_a, max(na_a)] real<lower=0> se;  // standard errors
  array[ns_a, max(na_a)] int<lower=1, upper=nt> t_a;
  
  real prior_mu_mean;
  real<lower=0> prior_mu_sd;
  real prior_d_mean;
  real<lower=0> prior_d_sd;
  
  int<lower=0, upper=1> has_metareg;
  array[ns_a, max(na_a)] real x_a;
  int<lower=0, upper=3> metareg_type;
  real prior_beta_mean;
  real<lower=0> prior_beta_sd;
}

transformed data {
  int max_arms = max(na_a);
}

parameters {
  vector[ns_a] mu;
  vector[nt-1] d_raw;
  vector[nt-1] beta_raw;
  real b;
  real<lower=0> gamma;
  real B;
}

transformed parameters {
  vector[nt] d;
  vector[nt] beta;
  
  d[1] = 0;
  d[2:nt] = d_raw;
  
  beta[1] = 0;
  if (has_metareg == 0) {
    beta[2:nt] = rep_vector(0, nt-1);
  } else if (metareg_type == 1) {
    beta[2:nt] = beta_raw;
  } else if (metareg_type == 2) {
    beta[2:nt] = b + gamma * beta_raw;
  } else if (metareg_type == 3) {
    beta[2:nt] = rep_vector(B, nt-1);
  }
}

model {
  mu ~ normal(prior_mu_mean, prior_mu_sd);
  d_raw ~ normal(prior_d_mean, prior_d_sd);
  
  if (has_metareg == 1) {
    if (metareg_type == 1) {
      beta_raw ~ student_t(1, prior_beta_mean, prior_beta_sd);
    } else if (metareg_type == 2) {
      beta_raw ~ std_normal();
      b ~ student_t(1, prior_beta_mean, prior_beta_sd);
      gamma ~ uniform(0, prior_beta_sd);
    } else if (metareg_type == 3) {
      B ~ student_t(1, prior_beta_mean, prior_beta_sd);
    }
  }
  
  for (i in 1:ns_a) {
    for (k in 1:na_a[i]) {
      real theta = mu[i] + d[t_a[i,k]] - d[t_a[i,1]];
      if (has_metareg == 1) {
        theta += (beta[t_a[i,k]] - beta[t_a[i,1]]) * x_a[i,k];
      }
      y[i,k] ~ normal(theta, se[i,k]);
    }
  }
}

generated quantities {
  array[ns_a, max_arms] real theta_a;
  array[ns_a, max_arms] real dev_a;
  real totresdev = 0;
  
  for (i in 1:ns_a) {
    for (k in 1:na_a[i]) {
      theta_a[i,k] = mu[i] + d[t_a[i,k]] - d[t_a[i,1]];
      if (has_metareg == 1) {
        theta_a[i,k] += (beta[t_a[i,k]] - beta[t_a[i,1]]) * x_a[i,k];
      }
      dev_a[i,k] = (y[i,k] - theta_a[i,k])^2 / (se[i,k]^2);
      totresdev += dev_a[i,k];
    }
  }
}

