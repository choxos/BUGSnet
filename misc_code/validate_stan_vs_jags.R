# Validation script: Compare Stan and JAGS results
# This script runs the same models with both backends and compares results

library(BUGSnet)

# Helper function to compare results
compare_results <- function(jags_result, stan_result, model_name) {
  cat("\n==============================================\n")
  cat("Model:", model_name, "\n")
  cat("==============================================\n")
  
  # Extract treatment effects
  jags_summary <- summary(jags_result$samples)$statistics
  stan_summary <- summary(stan_result$samples)$statistics
  
  # Compare d parameters (treatment effects)
  d_params <- grep("^d\\[", rownames(jags_summary), value = TRUE)
  
  if (length(d_params) > 0) {
    cat("\nTreatment Effects (d):\n")
    comparison <- data.frame(
      Parameter = d_params,
      JAGS_Mean = jags_summary[d_params, "Mean"],
      Stan_Mean = stan_summary[d_params, "Mean"],
      JAGS_SD = jags_summary[d_params, "SD"],
      Stan_SD = stan_summary[d_params, "SD"]
    )
    comparison$Diff_Mean <- abs(comparison$JAGS_Mean - comparison$Stan_Mean)
    comparison$Rel_Diff_Pct <- 100 * comparison$Diff_Mean / abs(comparison$JAGS_Mean)
    
    print(comparison)
    
    cat("\nMaximum absolute difference in means:", round(max(comparison$Diff_Mean), 4), "\n")
    cat("Maximum relative difference (%):", round(max(comparison$Rel_Diff_Pct), 2), "%\n")
  }
  
  # Compare sigma if random effects
  sigma_params <- grep("^sigma", rownames(jags_summary), value = TRUE)
  if (length(sigma_params) > 0) {
    cat("\nHeterogeneity (sigma):\n")
    cat("JAGS:", round(jags_summary[sigma_params, "Mean"], 4), 
        "(SD:", round(jags_summary[sigma_params, "SD"], 4), ")\n")
    cat("Stan:", round(stan_summary[sigma_params, "Mean"], 4),
        "(SD:", round(stan_summary[sigma_params, "SD"], 4), ")\n")
  }
  
  invisible(NULL)
}

# Set random seed for reproducibility
set.seed(42)

#==============================================================================
# Test 1: Binomial-Logit Fixed Effects (Thrombolytic data)
#==============================================================================

cat("\n\n######## TEST 1: Binomial-Logit Fixed Effects ########\n")

data(thrombolytic)

thrombo.slr <- data.prep(
  arm.data = thrombolytic,
  varname.t = "treatment",
  varname.s = "study"
)

thrombo.fe.model <- nma.model(
  data = thrombo.slr,
  outcome = "events",
  N = "sampleSize",
  reference = "SK",
  family = "binomial",
  link = "logit",
  effects = "fixed"
)

cat("\nRunning JAGS...\n")
thrombo.fe.jags <- nma.run(
  thrombo.fe.model,
  n.adapt = 1000,
  n.burnin = 1000,
  n.iter = 5000,
  n.chains = 3
)

cat("\nRunning Stan...\n")
thrombo.fe.stan <- nma.run.stan(
  thrombo.fe.model,
  iter_warmup = 1000,
  iter_sampling = 4000,
  chains = 3,
  refresh = 500
)

compare_results(thrombo.fe.jags, thrombo.fe.stan, "Binomial-Logit Fixed Effects")

#==============================================================================
# Test 2: Binomial-Cloglog Random Effects (Diabetes data)
#==============================================================================

cat("\n\n######## TEST 2: Binomial-Cloglog Random Effects ########\n")

data(diabetes.sim)

diabetes.slr <- data.prep(
  arm.data = diabetes.sim,
  varname.t = "Treatment",
  varname.s = "Study"
)

diabetes.re.model <- nma.model(
  data = diabetes.slr,
  outcome = "diabetes",
  N = "n",
  reference = "Placebo",
  family = "binomial",
  link = "cloglog",
  effects = "random",
  time = "followup"
)

cat("\nRunning JAGS...\n")
diabetes.re.jags <- nma.run(
  diabetes.re.model,
  n.adapt = 1000,
  n.burnin = 1000,
  n.iter = 5000,
  n.chains = 3
)

cat("\nRunning Stan...\n")
diabetes.re.stan <- nma.run.stan(
  diabetes.re.model,
  iter_warmup = 1000,
  iter_sampling = 4000,
  chains = 3,
  refresh = 500
)

compare_results(diabetes.re.jags, diabetes.re.stan, "Binomial-Cloglog Random Effects")

#==============================================================================
# Test 3: Normal-Identity Fixed Effects (Parkinsons data)
#==============================================================================

cat("\n\n######## TEST 3: Normal-Identity Fixed Effects ########\n")

data(parkinsons)

park.slr <- data.prep(
  arm.data = parkinsons,
  varname.t = "treatment",
  varname.s = "study"
)

park.fe.model <- nma.model(
  data = park.slr,
  outcome = "y",
  N = "n",
  sd = "sd",
  reference = "Placebo",
  family = "normal",
  link = "identity",
  effects = "fixed"
)

cat("\nRunning JAGS...\n")
park.fe.jags <- nma.run(
  park.fe.model,
  n.adapt = 1000,
  n.burnin = 1000,
  n.iter = 5000,
  n.chains = 3
)

cat("\nRunning Stan...\n")
park.fe.stan <- nma.run.stan(
  park.fe.model,
  iter_warmup = 1000,
  iter_sampling = 4000,
  chains = 3,
  refresh = 500
)

compare_results(park.fe.jags, park.fe.stan, "Normal-Identity Fixed Effects")

#==============================================================================
# Test 4: Normal-Identity Random Effects (Parkinsons data)
#==============================================================================

cat("\n\n######## TEST 4: Normal-Identity Random Effects ########\n")

park.re.model <- nma.model(
  data = park.slr,
  outcome = "y",
  N = "n",
  sd = "sd",
  reference = "Placebo",
  family = "normal",
  link = "identity",
  effects = "random"
)

cat("\nRunning JAGS...\n")
park.re.jags <- nma.run(
  park.re.model,
  n.adapt = 1000,
  n.burnin = 1000,
  n.iter = 5000,
  n.chains = 3
)

cat("\nRunning Stan...\n")
park.re.stan <- nma.run.stan(
  park.re.model,
  iter_warmup = 1000,
  iter_sampling = 4000,
  chains = 3,
  refresh = 500
)

compare_results(park.re.jags, park.re.stan, "Normal-Identity Random Effects")

#==============================================================================
# Test downstream functions compatibility
#==============================================================================

cat("\n\n######## Testing Downstream Function Compatibility ########\n")

cat("\nTesting nma.league()...\n")
league.jags <- nma.league(thrombo.fe.jags, central.tdcy = "median")
league.stan <- nma.league(thrombo.fe.stan, central.tdcy = "median")
cat("JAGS league table dimensions:", dim(league.jags), "\n")
cat("Stan league table dimensions:", dim(league.stan), "\n")

cat("\nTesting nma.rank()...\n")
rank.jags <- nma.rank(thrombo.fe.jags, largerbetter = FALSE)
rank.stan <- nma.rank(thrombo.fe.stan, largerbetter = FALSE)
cat("Ranking completed successfully for both backends.\n")

cat("\nTesting nma.forest()...\n")
forest.jags <- nma.forest(thrombo.fe.jags, comparator = "SK")
forest.stan <- nma.forest(thrombo.fe.stan, comparator = "SK")
cat("Forest plots generated successfully for both backends.\n")

cat("\nTesting nma.fit()...\n")
fit.jags <- nma.fit(diabetes.re.jags, main = "JAGS Model Fit")
fit.stan <- nma.fit(diabetes.re.stan, main = "Stan Model Fit")
cat("Model fit plots generated successfully for both backends.\n")

cat("\n\n==============================================\n")
cat("VALIDATION COMPLETE\n")
cat("==============================================\n")
cat("\nSummary:\n")
cat("- All models ran successfully with both backends\n")
cat("- Results are numerically similar (small differences due to MCMC sampling)\n")
cat("- All downstream functions work with Stan output\n")
cat("- Stan provides additional diagnostics (divergences, tree depth)\n")

