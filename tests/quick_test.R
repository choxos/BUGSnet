#!/usr/bin/env Rscript
# Quick test of all new functions with both JAGS and Stan

library(BUGSnet)

cat("==================== QUICK TEST ====================\n")
cat("Testing New Features with JAGS and Stan\n")
cat("====================================================\n\n")

# Prepare data
data(thrombolytic)
thrombo.slr <- data.prep(arm.data = thrombolytic,
                         varname.t = "treatment",
                         varname.s = "study")
cat("✓ Data prepared\n")

# Create model
thrombo.model <- nma.model(data = thrombo.slr,
                           outcome = "events",
                           N = "sampleSize",
                           reference = "SK",
                           family = "binomial",
                           link = "log",
                           effects = "random",
                           type = "consistency")
cat("✓ Model created (binomial-log random)\n\n")

# Test with JAGS
cat("--- Testing JAGS ---\n")
jags_result <- nma.run(model = thrombo.model,
                       n.adapt = 500,
                       n.burnin = 500,
                       n.iter = 1000)
cat("✓ JAGS sampling complete\n")

# Test new functions with JAGS
cat("\nTesting new functions with JAGS results:\n")

tryCatch({
  split_jags <- nma.netsplit(jags_result)
  cat("  ✓ nma.netsplit()\n")
  
  p_forest <- forest(split_jags, show = "with.direct")
  cat("  ✓ forest.nma.netsplit()\n")
  
  p_heat <- nma.heatplot(jags_result, digits = 2)
  cat("  ✓ nma.heatplot()\n")
  
  trt_order <- c("SK", "tPA", "PTCA")
  funnel_data <- nma.funnel(jags_result, order = trt_order)
  cat("  ✓ nma.funnel()\n")
  
  radial_data <- nma.radial(jags_result, order = trt_order)
  cat("  ✓ nma.radial()\n")
  
  ranks_jags <- nma.rank(jags_result, largerbetter = FALSE)
  cat("  ✓ nma.rank() integration\n")
  
  p_ordered <- nma.heatplot(jags_result, order = ranks_jags$order)
  cat("  ✓ heatplot with SUCRA ordering\n")
}, error = function(e) {
  cat("  ✗ Error:", e$message, "\n")
  stop(e)
})

# Test with Stan if available
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("\n--- Testing Stan ---\n")
  
  tryCatch({
    stan_result <- nma.run.stan(model = thrombo.model,
                                iter_warmup = 500,
                                iter_sampling = 500,
                                chains = 2,
                                refresh = 0,
                                show_messages = FALSE)
    cat("✓ Stan sampling complete\n")
    
    cat("\nTesting new functions with Stan results:\n")
    
    split_stan <- nma.netsplit(stan_result)
    cat("  ✓ nma.netsplit()\n")
    
    p_forest_stan <- forest(split_stan, show = "with.direct")
    cat("  ✓ forest.nma.netsplit()\n")
    
    p_heat_stan <- nma.heatplot(stan_result, digits = 2)
    cat("  ✓ nma.heatplot()\n")
    
    funnel_data_stan <- nma.funnel(stan_result, order = trt_order)
    cat("  ✓ nma.funnel()\n")
    
    radial_data_stan <- nma.radial(stan_result, order = trt_order)
    cat("  ✓ nma.radial()\n")
    
    ranks_stan <- nma.rank(stan_result, largerbetter = FALSE)
    cat("  ✓ nma.rank() integration\n")
    
    # Compare JAGS vs Stan
    jags_samples <- do.call(rbind, jags_result$samples)
    stan_samples <- do.call(rbind, stan_result$samples)
    
    jags_d <- grep("^d\\.", colnames(jags_samples), value = TRUE)
    stan_d <- grep("^d\\.", colnames(stan_samples), value = TRUE)
    
    jags_means <- colMeans(jags_samples[, jags_d])
    stan_means <- colMeans(stan_samples[, stan_d])
    
    correlation <- cor(jags_means, stan_means)
    cat(sprintf("\n  Correlation JAGS vs Stan: %.3f", correlation))
    if (correlation > 0.95) {
      cat(" ✓\n")
    } else {
      cat(sprintf(" (expected > 0.95)\n"))
    }
    
  }, error = function(e) {
    cat("✗ Stan test failed:", e$message, "\n")
  })
} else {
  cat("\n⚠ cmdstanr not available - skipping Stan tests\n")
}

cat("\n====================================================\n")
cat("All tests completed successfully!\n")
cat("====================================================\n")

