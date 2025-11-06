# Test Stan integration with BUGSnet

test_that("cmdstanr package check works", {
  # This should either succeed or give informative error
  expect_error(
    check_cmdstan_installation(),
    regexp = NA  # No error if cmdstanr is installed
  )
})

test_that("Stan model file selection works", {
  # Test that correct Stan files are selected
  expect_true(
    file.exists(select_stan_model("binomial", "logit", "fixed", "consistency"))
  )
  
  expect_true(
    file.exists(select_stan_model("binomial", "logit", "random", "consistency"))
  )
  
  expect_true(
    file.exists(select_stan_model("normal", "identity", "fixed", "consistency"))
  )
  
  expect_true(
    file.exists(select_stan_model("poisson", "log", "random", "consistency"))
  )
  
  # Test that invalid combinations give error
  expect_error(
    select_stan_model("invalid", "link", "fixed", "consistency"),
    "Stan model file not found"
  )
})

test_that("Data conversion works for binomial models", {
  skip_if_not_installed("cmdstanr")
  
  data(thrombolytic)
  
  thrombo.slr <- data.prep(
    arm.data = thrombolytic,
    varname.t = "treatment",
    varname.s = "study"
  )
  
  thrombo.model <- nma.model(
    data = thrombo.slr,
    outcome = "events",
    N = "sampleSize",
    reference = "SK",
    family = "binomial",
    link = "logit",
    effects = "fixed"
  )
  
  stan_data <- convert_bugs_to_stan_data(thrombo.model)
  
  # Check required fields
  expect_true("ns_a" %in% names(stan_data))
  expect_true("nt" %in% names(stan_data))
  expect_true("na_a" %in% names(stan_data))
  expect_true("r" %in% names(stan_data))
  expect_true("n" %in% names(stan_data))
  expect_true("t_a" %in% names(stan_data))
  
  # Check prior fields
  expect_true("prior_mu_mean" %in% names(stan_data))
  expect_true("prior_mu_sd" %in% names(stan_data))
  expect_true("prior_d_mean" %in% names(stan_data))
  expect_true("prior_d_sd" %in% names(stan_data))
})

test_that("nma.run.stan works with binomial-logit fixed effects", {
  skip_if_not_installed("cmdstanr")
  skip_if_not(cmdstanr::cmdstan_version(error_on_NA = FALSE) != "NA", 
              "CmdStan not installed")
  
  data(thrombolytic)
  
  thrombo.slr <- data.prep(
    arm.data = thrombolytic,
    varname.t = "treatment",
    varname.s = "study"
  )
  
  thrombo.model <- nma.model(
    data = thrombo.slr,
    outcome = "events",
    N = "sampleSize",
    reference = "SK",
    family = "binomial",
    link = "logit",
    effects = "fixed"
  )
  
  # Run with minimal iterations for testing
  result <- nma.run.stan(
    thrombo.model,
    iter_warmup = 200,
    iter_sampling = 200,
    chains = 2,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Check result structure
  expect_s3_class(result, "BUGSnetRun")
  expect_true("samples" %in% names(result))
  expect_true("model" %in% names(result))
  expect_true("stan_fit" %in% names(result))
  
  # Check samples are in coda format
  expect_s3_class(result$samples, "mcmc.list")
})

test_that("Downstream functions work with Stan results", {
  skip_if_not_installed("cmdstanr")
  skip_if_not(cmdstanr::cmdstan_version(error_on_NA = FALSE) != "NA", 
              "CmdStan not installed")
  
  data(thrombolytic)
  
  thrombo.slr <- data.prep(
    arm.data = thrombolytic,
    varname.t = "treatment",
    varname.s = "study"
  )
  
  thrombo.model <- nma.model(
    data = thrombo.slr,
    outcome = "events",
    N = "sampleSize",
    reference = "SK",
    family = "binomial",
    link = "logit",
    effects = "fixed"
  )
  
  result <- nma.run.stan(
    thrombo.model,
    iter_warmup = 200,
    iter_sampling = 200,
    chains = 2,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Test nma.fit
  expect_error(
    fit_result <- nma.fit(result, main = "Test"),
    NA
  )
  
  # Test nma.league
  expect_error(
    league_result <- nma.league(result, central.tdcy = "median"),
    NA
  )
  expect_true(is.data.frame(league_result) || is.matrix(league_result))
  
  # Test nma.rank
  expect_error(
    rank_result <- nma.rank(result, largerbetter = FALSE),
    NA
  )
  
  # Test nma.forest
  expect_error(
    forest_result <- nma.forest(result, comparator = "SK"),
    NA
  )
})

test_that("Stan and JAGS give similar results", {
  skip_if_not_installed("cmdstanr")
  skip_if_not(cmdstanr::cmdstan_version(error_on_NA = FALSE) != "NA", 
              "CmdStan not installed")
  
  data(thrombolytic)
  
  thrombo.slr <- data.prep(
    arm.data = thrombolytic,
    varname.t = "treatment",
    varname.s = "study"
  )
  
  thrombo.model <- nma.model(
    data = thrombo.slr,
    outcome = "events",
    N = "sampleSize",
    reference = "SK",
    family = "binomial",
    link = "logit",
    effects = "fixed"
  )
  
  # Run with JAGS
  set.seed(123)
  result_jags <- nma.run(
    thrombo.model,
    n.adapt = 100,
    n.burnin = 200,
    n.iter = 500,
    n.chains = 2
  )
  
  # Run with Stan
  set.seed(123)
  result_stan <- nma.run.stan(
    thrombo.model,
    iter_warmup = 200,
    iter_sampling = 300,
    chains = 2,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Extract treatment effects
  d_jags <- summary(result_jags$samples)$statistics[grep("^d\\[", rownames(summary(result_jags$samples)$statistics)), "Mean"]
  d_stan <- summary(result_stan$samples)$statistics[grep("^d\\[", rownames(summary(result_stan$samples)$statistics)), "Mean"]
  
  # Results should be similar (within reasonable Monte Carlo error)
  # Using generous tolerance for test stability
  expect_equal(d_jags, d_stan, tolerance = 0.5)
})

test_that("Random effects models work with Stan", {
  skip_if_not_installed("cmdstanr")
  skip_if_not(cmdstanr::cmdstan_version(error_on_NA = FALSE) != "NA", 
              "CmdStan not installed")
  
  data(diabetes.sim)
  
  diabetes.slr <- data.prep(
    arm.data = diabetes.sim,
    varname.t = "Treatment",
    varname.s = "Study"
  )
  
  diabetes.model <- nma.model(
    data = diabetes.slr,
    outcome = "diabetes",
    N = "n",
    reference = "Placebo",
    family = "binomial",
    link = "cloglog",
    effects = "random",
    time = "followup"
  )
  
  result <- nma.run.stan(
    diabetes.model,
    iter_warmup = 200,
    iter_sampling = 200,
    chains = 2,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Check that sigma is monitored
  samples_summary <- summary(result$samples)$statistics
  expect_true(any(grepl("^sigma", rownames(samples_summary))))
})

