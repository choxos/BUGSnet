#' Check if cmdstanr is installed and configured
#' @noRd
check_cmdstan_installation <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required but not installed.\n",
         "Install it with:\n",
         "  install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))\n",
         "After installing cmdstanr, you may need to install CmdStan with:\n",
         "  cmdstanr::install_cmdstan()",
         call. = FALSE)
  }
  
  # Check if CmdStan is installed
  cmdstan_ver <- tryCatch({
    cmdstanr::cmdstan_version(error_on_NA = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(cmdstan_ver) || is.na(cmdstan_ver)) {
    message("CmdStan not found. Attempting to install...")
    tryCatch({
      cmdstanr::install_cmdstan()
    }, error = function(e) {
      stop("Failed to install CmdStan: ", e$message, "\n",
           "You can try installing manually with: cmdstanr::install_cmdstan()",
           call. = FALSE)
    })
  }
  
  invisible(TRUE)
}


#' Select appropriate Stan model file based on model characteristics
#' @noRd
select_stan_model <- function(family, link, effects, type) {
  # Map family/link combinations to model names
  model_name <- paste(family, link, effects, type, sep = "_")
  
  # Construct file path
  stan_file <- system.file("stan", paste0(model_name, ".stan"), package = "BUGSnet")
  
  if (stan_file == "" || !file.exists(stan_file)) {
    stop("Stan model file not found for combination: ", model_name,
         "\nFamily: ", family, ", Link: ", link, 
         ", Effects: ", effects, ", Type: ", type,
         call. = FALSE)
  }
  
  return(stan_file)
}


#' Convert BUGSnetModel data to Stan data format
#' @noRd
convert_bugs_to_stan_data <- function(model) {
  bugsdata <- model$data
  family <- model$family
  link <- model$link
  effects <- model$effects
  
  # Calculate max number of arms
  max_arms <- max(bugsdata$na_a)
  
  # Initialize Stan data list
  stan_data <- list(
    ns_a = bugsdata$ns_a,
    nt = bugsdata$nt,
    na_a = bugsdata$na_a
  )
  
  # Add treatment indices - replace NAs with 1 (will be ignored by Stan based on na_a)
  stan_data$t_a <- bugsdata$t_a
  stan_data$t_a[is.na(stan_data$t_a)] <- 1
  
  # Add family-specific data - replace NAs with 0 (will be ignored by Stan based on na_a)
  if (family == "binomial") {
    stan_data$r <- bugsdata$r
    stan_data$r[is.na(stan_data$r)] <- 0
    stan_data$n <- bugsdata$n
    stan_data$n[is.na(stan_data$n)] <- 1  # avoid division by zero
    
    # Add time for cloglog link
    if (link == "cloglog") {
      stan_data$time <- bugsdata$time
      stan_data$time[is.na(stan_data$time)] <- 1
    }
  } else if (family == "normal") {
    stan_data$y <- bugsdata$y
    stan_data$y[is.na(stan_data$y)] <- 0
    stan_data$se <- bugsdata$se
    stan_data$se[is.na(stan_data$se)] <- 1  # avoid division by zero
  } else if (family == "poisson") {
    stan_data$r <- bugsdata$r
    stan_data$r[is.na(stan_data$r)] <- 0
    stan_data$E <- bugsdata$E
    stan_data$E[is.na(stan_data$E)] <- 1  # avoid division by zero
  }
  
  # Add priors - use reasonable defaults based on scale
  # These match what BUGSnet uses for default priors
  max.delta <- switch(model$scale,
                     "Odds Ratio" = 5,
                     "Risk Ratio" = 2,
                     "Mean Difference" = 10,
                     "Hazard Ratio" = 3,
                     "Rate Ratio" = 3,
                     5)  # default
  
  # Convert prior precision to SD (JAGS uses precision, Stan uses SD)
  if (model$prior.mu == "DEFAULT") {
    stan_data$prior_mu_mean <- 0
    stan_data$prior_mu_sd <- 15 * max.delta
  } else {
    # Parse custom prior
    stan_data$prior_mu_mean <- 0
    stan_data$prior_mu_sd <- 15 * max.delta
  }
  
  if (model$prior.d == "DEFAULT") {
    stan_data$prior_d_mean <- 0
    stan_data$prior_d_sd <- 15 * max.delta
  } else {
    stan_data$prior_d_mean <- 0
    stan_data$prior_d_sd <- 15 * max.delta
  }
  
  if (effects == "random") {
    if (model$prior.sigma == "DEFAULT") {
      stan_data$prior_sigma_max <- max.delta
    } else {
      stan_data$prior_sigma_max <- max.delta
    }
  }
  
  # Add meta-regression data
  if (!is.null(model$covariate)) {
    stan_data$has_metareg <- 1
    stan_data$x_a <- bugsdata$x_a
    stan_data$x_a[is.na(stan_data$x_a)] <- 0
    
    # Determine meta-regression type
    if (model$prior.beta == "UNRELATED") {
      stan_data$metareg_type <- 1
    } else if (model$prior.beta == "EXCHANGEABLE") {
      stan_data$metareg_type <- 2
    } else if (model$prior.beta == "EQUAL") {
      stan_data$metareg_type <- 3
    } else {
      stan_data$metareg_type <- 1  # default to UNRELATED
    }
    
    stan_data$prior_beta_mean <- 0
    stan_data$prior_beta_sd <- max.delta
  } else {
    stan_data$has_metareg <- 0
    # Still need to provide placeholder values for Stan
    stan_data$x_a <- matrix(0, nrow = stan_data$ns_a, ncol = max_arms)
    stan_data$metareg_type <- 0
    stan_data$prior_beta_mean <- 0
    stan_data$prior_beta_sd <- 1
  }
  
  return(stan_data)
}


#' Compile Stan model
#' @noRd
compile_stan_model <- function(stan_file) {
  check_cmdstan_installation()
  
  tryCatch({
    mod <- cmdstanr::cmdstan_model(stan_file)
    return(mod)
  }, error = function(e) {
    stop("Failed to compile Stan model: ", e$message, call. = FALSE)
  })
}


#' Convert cmdstanr fit to BUGSnetRun format
#' @importFrom coda as.mcmc as.mcmc.list
#' @noRd
convert_stan_to_bugsnet <- function(fit, model) {
  # Extract draws from cmdstanr fit
  draws_array <- fit$draws(format = "array")
  
  # Convert to coda mcmc.list format for compatibility
  n_chains <- dim(draws_array)[2]
  n_iter <- dim(draws_array)[1]
  
  # Get parameter names
  param_names <- dimnames(draws_array)[[3]]
  
  # Create mcmc.list
  mcmc_list <- list()
  for (chain in 1:n_chains) {
    chain_draws <- draws_array[, chain, , drop = FALSE]
    dim(chain_draws) <- c(n_iter, length(param_names))
    colnames(chain_draws) <- param_names
    mcmc_list[[chain]] <- coda::as.mcmc(chain_draws)
  }
  mcmc_list <- coda::as.mcmc.list(mcmc_list)
  
  # Create BUGSnetRun object
  brun <- structure(list(
    samples = mcmc_list,
    model = model,
    scale = model$scale,
    family = model$family,
    link = model$link,
    trt.key = as.character(t(model$trt.key[1])),
    stan_fit = fit  # Keep original Stan fit for diagnostics
  ), class = "BUGSnetRun")
  
  return(brun)
}


#' Get Stan sampling diagnostics
#' @noRd
get_stan_diagnostics <- function(fit) {
  diag <- fit$diagnostic_summary()
  
  issues <- list()
  
  # Check for divergences
  if (!is.null(diag$num_divergent) && sum(diag$num_divergent) > 0) {
    issues$divergences <- sum(diag$num_divergent)
  }
  
  # Check for max treedepth
  if (!is.null(diag$num_max_treedepth) && sum(diag$num_max_treedepth) > 0) {
    issues$max_treedepth <- sum(diag$num_max_treedepth)
  }
  
  # Get summary for Rhat and ESS
  summary_fit <- fit$summary()
  
  # Check Rhat
  rhat_vals <- summary_fit$rhat
  max_rhat <- max(rhat_vals, na.rm = TRUE)
  if (max_rhat > 1.05) {
    issues$max_rhat <- max_rhat
  }
  
  # Check ESS
  ess_bulk <- summary_fit$ess_bulk
  min_ess <- min(ess_bulk, na.rm = TRUE)
  if (min_ess < 400) {
    issues$min_ess <- min_ess
  }
  
  return(issues)
}


#' Print diagnostic warnings
#' @noRd
print_stan_warnings <- function(diagnostics, sampling_args) {
  if (length(diagnostics) == 0) {
    message("All MCMC diagnostics look good.")
    return(invisible(NULL))
  }
  
  if (!is.null(diagnostics$divergences)) {
    warning(diagnostics$divergences, " divergent transitions detected. ",
            "Consider increasing adapt_delta (currently ", 
            sampling_args$adapt_delta, ")", call. = FALSE)
  }
  
  if (!is.null(diagnostics$max_treedepth)) {
    warning(diagnostics$max_treedepth, " iterations hit max treedepth. ",
            "Consider increasing max_treedepth (currently ", 
            sampling_args$max_treedepth, ")", call. = FALSE)
  }
  
  if (!is.null(diagnostics$max_rhat)) {
    warning("Some Rhat values > 1.05 (max = ", round(diagnostics$max_rhat, 3), "). ",
            "Chains may not have converged. Consider running more iterations.",
            call. = FALSE)
  }
  
  if (!is.null(diagnostics$min_ess)) {
    warning("Some bulk ESS values < 400 (min = ", round(diagnostics$min_ess, 1), "). ",
            "Consider running more iterations.", call. = FALSE)
  }
  
  invisible(NULL)
}

