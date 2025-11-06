# Comprehensive Test Script for netmeta-style Features
# Tests all new plotting functions with both JAGS and Stan backends

library(BUGSnet)

cat("=================================================================\n")
cat("Testing netmeta-style features in BUGSnet\n")
cat("Testing with both JAGS and Stan backends\n")
cat("=================================================================\n\n")

# Helper function to print test results
test_result <- function(test_name, result) {
  if (result) {
    cat(sprintf("✓ PASS: %s\n", test_name))
  } else {
    cat(sprintf("✗ FAIL: %s\n", test_name))
  }
}

# =================================================================
# Test 1: Data Preparation
# =================================================================
cat("\n--- Test 1: Data Preparation ---\n")
data(thrombolytic)
tryCatch({
  thrombo.slr <- data.prep(arm.data = thrombolytic,
                           varname.t = "treatment",
                           varname.s = "study")
  test_result("Data preparation", TRUE)
}, error = function(e) {
  test_result("Data preparation", FALSE)
  cat("Error:", e$message, "\n")
})

# =================================================================
# Test 2: Model Creation
# =================================================================
cat("\n--- Test 2: Model Creation ---\n")
tryCatch({
  thrombo.model <- nma.model(data = thrombo.slr,
                             outcome = "events",
                             N = "sampleSize",
                             reference = "SK",
                             family = "binomial",
                             link = "log",
                             effects = "random",
                             type = "consistency")
  test_result("Model creation", TRUE)
}, error = function(e) {
  test_result("Model creation", FALSE)
  cat("Error:", e$message, "\n")
})

# =================================================================
# Test 3: JAGS Backend
# =================================================================
cat("\n--- Test 3: JAGS Backend ---\n")
cat("Running MCMC with JAGS (this may take a minute)...\n")
jags_result <- NULL
tryCatch({
  jags_result <- nma.run(model = thrombo.model,
                         n.adapt = 500,
                         n.burnin = 500,
                         n.iter = 1000)
  test_result("JAGS sampling", !is.null(jags_result))
}, error = function(e) {
  test_result("JAGS sampling", FALSE)
  cat("Error:", e$message, "\n")
})

if (!is.null(jags_result)) {
  # Test new functions with JAGS results
  
  # Test 3.1: nma.netsplit
  cat("\n  Testing nma.netsplit() with JAGS...\n")
  tryCatch({
    split_jags <- nma.netsplit(jags_result)
    test_result("  nma.netsplit() with JAGS", !is.null(split_jags))
    
    # Check structure
    has_comparison <- "comparison" %in% names(split_jags)
    has_nma <- "nma" %in% names(split_jags)
    has_direct <- "direct" %in% names(split_jags)
    has_indirect <- "indirect" %in% names(split_jags)
    test_result("  netsplit structure correct", 
                has_comparison && has_nma && has_direct && has_indirect)
    
  }, error = function(e) {
    test_result("  nma.netsplit() with JAGS", FALSE)
    cat("  Error:", e$message, "\n")
  })
  
  # Test 3.2: forest.nma.netsplit
  cat("\n  Testing forest.nma.netsplit() with JAGS...\n")
  tryCatch({
    if (exists("split_jags")) {
      p_forest <- forest(split_jags, show = "with.direct")
      test_result("  forest.nma.netsplit() with JAGS", !is.null(p_forest))
      test_result("  forest plot is ggplot2", inherits(p_forest, "gg"))
    }
  }, error = function(e) {
    test_result("  forest.nma.netsplit() with JAGS", FALSE)
    cat("  Error:", e$message, "\n")
  })
  
  # Test 3.3: nma.heatplot
  cat("\n  Testing nma.heatplot() with JAGS...\n")
  tryCatch({
    p_heat <- nma.heatplot(jags_result, digits = 2)
    test_result("  nma.heatplot() with JAGS", !is.null(p_heat))
    test_result("  heatplot is ggplot2", inherits(p_heat, "gg"))
  }, error = function(e) {
    test_result("  nma.heatplot() with JAGS", FALSE)
    cat("  Error:", e$message, "\n")
  })
  
  # Test 3.4: nma.funnel
  cat("\n  Testing nma.funnel() with JAGS...\n")
  tryCatch({
    trt_order <- c("SK", "tPA", "PTCA")
    funnel_data <- nma.funnel(jags_result, order = trt_order)
    test_result("  nma.funnel() with JAGS", !is.null(funnel_data))
  }, error = function(e) {
    test_result("  nma.funnel() with JAGS", FALSE)
    cat("  Error:", e$message, "\n")
  })
  
  # Test 3.5: nma.radial
  cat("\n  Testing nma.radial() with JAGS...\n")
  tryCatch({
    trt_order <- c("SK", "tPA", "PTCA")
    radial_data <- nma.radial(jags_result, order = trt_order)
    test_result("  nma.radial() with JAGS", !is.null(radial_data))
  }, error = function(e) {
    test_result("  nma.radial() with JAGS", FALSE)
    cat("  Error:", e$message, "\n")
  })
  
  # Test 3.6: Integration with existing nma.rank
  cat("\n  Testing integration with nma.rank()...\n")
  tryCatch({
    ranks_jags <- nma.rank(jags_result, largerbetter = FALSE)
    test_result("  nma.rank() with JAGS", !is.null(ranks_jags))
    
    has_ranktable <- "ranktable" %in% names(ranks_jags)
    has_sucra <- "sucratable" %in% names(ranks_jags)
    has_plots <- "sucraplot" %in% names(ranks_jags) && "rankogram" %in% names(ranks_jags)
    test_result("  rank structure correct", 
                has_ranktable && has_sucra && has_plots)
    
    # Test using rank order in heatplot
    if ("order" %in% names(ranks_jags)) {
      p_ordered <- nma.heatplot(jags_result, order = ranks_jags$order)
      test_result("  heatplot with SUCRA ordering", !is.null(p_ordered))
    }
  }, error = function(e) {
    test_result("  nma.rank() with JAGS", FALSE)
    cat("  Error:", e$message, "\n")
  })
  
  # Test 3.7: Existing functions still work
  cat("\n  Testing backward compatibility with existing functions...\n")
  tryCatch({
    league_jags <- nma.league(jags_result)
    test_result("  nma.league() still works", !is.null(league_jags))
    
    # Wrap forest plot in pdf to avoid display issues
    pdf(NULL)
    forest_test <- nma.forest(jags_result, comparator = "SK")
    dev.off()
    test_result("  nma.forest() still works", TRUE)
    
  }, error = function(e) {
    test_result("  Backward compatibility", FALSE)
    cat("  Error:", e$message, "\n")
  })
}

# =================================================================
# Test 4: Stan Backend (if available)
# =================================================================
cat("\n--- Test 4: Stan Backend (cmdstanr) ---\n")

stan_available <- FALSE
tryCatch({
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    stan_available <- TRUE
    cat("cmdstanr is installed\n")
  } else {
    cat("cmdstanr is not installed - skipping Stan tests\n")
  }
}, error = function(e) {
  cat("cmdstanr is not available - skipping Stan tests\n")
})

if (stan_available) {
  cat("Running MCMC with Stan (this may take a minute)...\n")
  stan_result <- NULL
  tryCatch({
    stan_result <- nma.run.stan(model = thrombo.model,
                                iter_warmup = 500,
                                iter_sampling = 500,
                                chains = 2,
                                refresh = 0,
                                show_messages = FALSE)
    test_result("Stan sampling", !is.null(stan_result))
  }, error = function(e) {
    test_result("Stan sampling", FALSE)
    cat("Error:", e$message, "\n")
    cat("Note: This may fail if CmdStan is not properly installed\n")
  })
  
  if (!is.null(stan_result)) {
    # Test new functions with Stan results
    
    # Test 4.1: nma.netsplit
    cat("\n  Testing nma.netsplit() with Stan...\n")
    tryCatch({
      split_stan <- nma.netsplit(stan_result)
      test_result("  nma.netsplit() with Stan", !is.null(split_stan))
    }, error = function(e) {
      test_result("  nma.netsplit() with Stan", FALSE)
      cat("  Error:", e$message, "\n")
    })
    
    # Test 4.2: forest.nma.netsplit
    cat("\n  Testing forest.nma.netsplit() with Stan...\n")
    tryCatch({
      if (exists("split_stan")) {
        p_forest_stan <- forest(split_stan, show = "with.direct")
        test_result("  forest.nma.netsplit() with Stan", !is.null(p_forest_stan))
      }
    }, error = function(e) {
      test_result("  forest.nma.netsplit() with Stan", FALSE)
      cat("  Error:", e$message, "\n")
    })
    
    # Test 4.3: nma.heatplot
    cat("\n  Testing nma.heatplot() with Stan...\n")
    tryCatch({
      p_heat_stan <- nma.heatplot(stan_result, digits = 2)
      test_result("  nma.heatplot() with Stan", !is.null(p_heat_stan))
    }, error = function(e) {
      test_result("  nma.heatplot() with Stan", FALSE)
      cat("  Error:", e$message, "\n")
    })
    
    # Test 4.4: nma.funnel
    cat("\n  Testing nma.funnel() with Stan...\n")
    tryCatch({
      trt_order <- c("SK", "tPA", "PTCA")
      funnel_data_stan <- nma.funnel(stan_result, order = trt_order)
      test_result("  nma.funnel() with Stan", !is.null(funnel_data_stan))
    }, error = function(e) {
      test_result("  nma.funnel() with Stan", FALSE)
      cat("  Error:", e$message, "\n")
    })
    
    # Test 4.5: nma.radial
    cat("\n  Testing nma.radial() with Stan...\n")
    tryCatch({
      trt_order <- c("SK", "tPA", "PTCA")
      radial_data_stan <- nma.radial(stan_result, order = trt_order)
      test_result("  nma.radial() with Stan", !is.null(radial_data_stan))
    }, error = function(e) {
      test_result("  nma.radial() with Stan", FALSE)
      cat("  Error:", e$message, "\n")
    })
    
    # Test 4.6: Integration with nma.rank
    cat("\n  Testing integration with nma.rank()...\n")
    tryCatch({
      ranks_stan <- nma.rank(stan_result, largerbetter = FALSE)
      test_result("  nma.rank() with Stan", !is.null(ranks_stan))
      
      # Test using rank order in heatplot
      if ("order" %in% names(ranks_stan)) {
        p_ordered_stan <- nma.heatplot(stan_result, order = ranks_stan$order)
        test_result("  heatplot with SUCRA ordering (Stan)", !is.null(p_ordered_stan))
      }
    }, error = function(e) {
      test_result("  nma.rank() with Stan", FALSE)
      cat("  Error:", e$message, "\n")
    })
    
    # Test 4.7: Existing functions still work with Stan
    cat("\n  Testing backward compatibility with Stan...\n")
    tryCatch({
      league_stan <- nma.league(stan_result)
      test_result("  nma.league() with Stan", !is.null(league_stan))
      
      pdf(NULL)
      forest_test_stan <- nma.forest(stan_result, comparator = "SK")
      dev.off()
      test_result("  nma.forest() with Stan", TRUE)
      
    }, error = function(e) {
      test_result("  Backward compatibility (Stan)", FALSE)
      cat("  Error:", e$message, "\n")
    })
    
    # Test 4.8: Compare JAGS vs Stan results
    if (!is.null(jags_result) && !is.null(stan_result)) {
      cat("\n  Comparing JAGS vs Stan results...\n")
      tryCatch({
        # Extract d parameters from both
        jags_samples <- do.call(rbind, jags_result$samples)
        stan_samples <- do.call(rbind, stan_result$samples)
        
        # Check that same parameters exist
        jags_d <- grep("^d\\.", colnames(jags_samples), value = TRUE)
        stan_d <- grep("^d\\.", colnames(stan_samples), value = TRUE)
        
        test_result("  Same parameters in JAGS and Stan", 
                    length(jags_d) == length(stan_d))
        
        # Check that estimates are in similar range (rough check)
        jags_means <- colMeans(jags_samples[, jags_d])
        stan_means <- colMeans(stan_samples[, stan_d])
        
        correlation <- cor(jags_means, stan_means)
        test_result("  JAGS and Stan results correlated", correlation > 0.95)
        
        cat(sprintf("    Correlation between JAGS and Stan: %.3f\n", correlation))
        
      }, error = function(e) {
        test_result("  JAGS vs Stan comparison", FALSE)
        cat("  Error:", e$message, "\n")
      })
    }
  }
}

# =================================================================
# Test 5: Error Handling
# =================================================================
cat("\n--- Test 5: Error Handling ---\n")

# Test with invalid inputs
cat("  Testing error handling...\n")

tryCatch({
  # Should fail gracefully
  nma.netsplit("not_a_BUGSnetRun_object")
  test_result("  Invalid input handling", FALSE)
}, error = function(e) {
  test_result("  Invalid input handling", TRUE)
})

if (!is.null(jags_result)) {
  tryCatch({
    # Should fail gracefully - order is mandatory
    nma.funnel(jags_result)
    test_result("  Missing mandatory argument handling", FALSE)
  }, error = function(e) {
    test_result("  Missing mandatory argument handling", TRUE)
  })
}

# =================================================================
# Test 6: Print and Summary Methods
# =================================================================
cat("\n--- Test 6: Print and Summary Methods ---\n")

if (!is.null(jags_result) && exists("split_jags")) {
  tryCatch({
    # Capture print output
    capture.output(print(split_jags))
    test_result("print.nma.netsplit()", TRUE)
  }, error = function(e) {
    test_result("print.nma.netsplit()", FALSE)
    cat("Error:", e$message, "\n")
  })
}

# =================================================================
# Summary
# =================================================================
cat("\n=================================================================\n")
cat("Test Suite Complete\n")
cat("=================================================================\n\n")

if (stan_available) {
  cat("✓ Tested with both JAGS and Stan backends\n")
} else {
  cat("⚠ Tested with JAGS only (cmdstanr not available)\n")
  cat("  To test Stan: install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))\n")
}

cat("\nAll new functions tested:\n")
cat("  • nma.netsplit() - Evidence splitting\n")
cat("  • forest.nma.netsplit() - Forest plots for split evidence\n")
cat("  • nma.heatplot() - Color-coded league tables\n")
cat("  • nma.funnel() - Comparison-adjusted funnel plots\n")
cat("  • nma.radial() - Comparison-adjusted radial plots\n")
cat("  • Integration with nma.rank() - Treatment rankings\n")
cat("  • Backward compatibility with existing functions\n")

cat("\n=================================================================\n")

