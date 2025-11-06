// Random effects inconsistency model for binomial data with log link

data {
  int<lower=0> ns_a;
  int<lower=0> nt;
  array[ns_a] int<lower=2> na_a;
  array[ns_a, max(na_a)] int<lower=0> r;
  array[ns_a, max(na_a)] int<lower=0> n;
  array[ns_a, max(na_a)] int<lower=1, upper=nt> t_a;
  
  real prior_d_mean;
  real<lower=0> prior_d_sd;
  real<lower=0> prior_sigma_max;
  
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
  vector<lower=0, upper=1>[ns_a] p_base;
  matrix[nt, nt] d_raw;
  real<lower=0> sigma;
  array[ns_a, max_arms] real delta_raw;
  vector[nt-1] beta_raw;
  real b;
  real<lower=0> gamma;
  real B;
}

transformed parameters {
  matrix[nt, nt] d;
  vector[nt] beta;
  vector[ns_a] mu = log(p_base);
  array[ns_a, max_arms] real delta;
  real sigma2 = sigma^2;
  
  for (i in 1:nt) d[i,i] = 0;
  for (c in 1:(nt-1)) {
    for (k in (c+1):nt) {
      d[c,k] = d_raw[c,k];
      d[k,c] = -d[c,k];
    }
  }
  
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
  
  for (i in 1:ns_a) {
    delta[i,1] = 0;
    for (k in 2:na_a[i]) {
      delta[i,k] = d[t_a[i,1], t_a[i,k]] + sigma * delta_raw[i,k];
    }
  }
}

model {
  p_base ~ uniform(0, 1);
  sigma ~ uniform(0, prior_sigma_max);
  
  for (c in 1:(nt-1)) {
    for (k in (c+1):nt) {
      d_raw[c,k] ~ normal(prior_d_mean, prior_d_sd);
    }
  }
  
  for (i in 1:ns_a) {
    for (k in 2:na_a[i]) {
      delta_raw[i,k] ~ std_normal();
    }
  }
  
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
      real log_p = mu[i] + delta[i,k];
      if (has_metareg == 1) {
        log_p += (beta[t_a[i,k]] - beta[t_a[i,1]]) * x_a[i,k];
      }
      // For log link: p = exp(log_p), constrained to [0,1]
      real p = fmin(exp(log_p), 0.9999);
      r[i,k] ~ binomial(n[i,k], p);
    }
  }
}

generated quantities {
  array[ns_a, max_arms] real rhat;
  array[ns_a, max_arms] real dev_a;
  real totresdev = 0;
  real tau = 1 / sigma2;
  
  for (i in 1:ns_a) {
    for (k in 1:na_a[i]) {
      real log_p = mu[i] + delta[i,k];
      real p;
      if (has_metareg == 1) {
        log_p += (beta[t_a[i,k]] - beta[t_a[i,1]]) * x_a[i,k];
      }
      p = exp(log_p);
      rhat[i,k] = p * n[i,k];
      
      if (r[i,k] > 0 && r[i,k] < n[i,k]) {
        dev_a[i,k] = 2 * (r[i,k] * (log(r[i,k]) - log(rhat[i,k])) + 
                          (n[i,k] - r[i,k]) * (log(n[i,k] - r[i,k]) - log(n[i,k] - rhat[i,k])));
      } else if (r[i,k] == 0) {
        dev_a[i,k] = 2 * (n[i,k] - r[i,k]) * (log(n[i,k] - r[i,k]) - log(n[i,k] - rhat[i,k]));
      } else {
        dev_a[i,k] = 2 * r[i,k] * (log(r[i,k]) - log(rhat[i,k]));
      }
      totresdev += dev_a[i,k];
    }
  }
}

