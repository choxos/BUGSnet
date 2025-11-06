#' @title
#' Run NMA model using Stan
#' 
#' @description
#' Takes a model from an object produced by \code{nma.model} and runs it using Stan via \code{cmdstanr}.
#' This is an alternative to \code{nma.run} which uses JAGS.
#' 
#' @param model A \code{BUGSnetModel} object produced by running \code{nma.model}.
#' @param monitor A vector of all variables that you would like to monitor. Default is "DEFAULT" which will monitor the relative treatment effects \code{d} 
#' as well as \code{sigma} when a random effects model is used and the regression coefficients \code{beta} when meta-regression is used.
#' @param iter_warmup Number of warmup (burn-in) iterations per chain. Default is 1000.
#' @param iter_sampling Number of post-warmup sampling iterations per chain. Default is 1000.
#' @param chains Number of MCMC chains. Default is 4.
#' @param parallel_chains Number of chains to run in parallel. Default is the same as \code{chains}.
#' @param thin Thinning rate. Default is 1 (no thinning).
#' @param adapt_delta Target average acceptance probability for adaptation. Default is 0.95. 
#' Increase if you see divergent transitions.
#' @param max_treedepth Maximum tree depth for the NUTS sampler. Default is 15. 
#' Increase if you see warnings about hitting max treedepth.
#' @param seed Random seed for reproducibility.
#' @param refresh How often to print progress updates. 0 suppresses output. Default is 100.
#' @param show_messages Logical. Whether to show compilation messages. Default is FALSE.
#' @param ... Additional arguments passed to \code{cmdstanr::sample()}.
#' 
#' @return \code{nma.run.stan} returns an object of class \code{BUGSnetRun} which is a list containing the following components:
#' @return \code{samples} - The MCMC samples in coda format, compatible with downstream BUGSnet functions.
#' @return \code{model} - The \code{BUGSnetModel} object obtained from \code{nma.model}.
#' @return \code{scale} - The scale of the outcome, based on the chosen family and link function.
#' @return \code{trt.key} - Treatments mapped to numbers, used in the model.
#' @return \code{family} - Family that was used for the NMA model (e.g normal, binomial, poisson)
#' @return \code{link} - Link function that was used for the NMA model (e.g identity, logit, log, cloglog)
#' @return \code{stan_fit} - The original cmdstanr fit object for advanced diagnostics.
#' 
#' @details
#' This function requires the \code{cmdstanr} package, which is not on CRAN. Install it with:
#' \code{install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))}
#' 
#' After installing cmdstanr, you may need to install CmdStan with \code{cmdstanr::install_cmdstan()}.
#' 
#' Stan uses Hamiltonian Monte Carlo (HMC) instead of Gibbs sampling, which can be more efficient
#' but requires compilation of the model. The first run may take longer due to compilation,
#' but subsequent runs with the same model type will be faster.
#' 
#' Key differences from \code{nma.run} (JAGS):
#' \itemize{
#'   \item Uses \code{iter_warmup} and \code{iter_sampling} instead of \code{n.adapt}, \code{n.burnin}, and \code{n.iter}
#'   \item Includes Stan-specific tuning parameters \code{adapt_delta} and \code{max_treedepth}
#'   \item Provides additional diagnostics like divergent transitions and tree depth warnings
#'   \item Generally faster sampling but requires initial compilation
#' }
#' 
#' @seealso
#' \code{\link{nma.model}}, \code{\link{nma.run}}, \code{\link{nma.fit}}, \code{\link{nma.league}}, 
#' \code{\link{nma.rank}}, \code{\link{nma.forest}}, \code{\link{nma.trace}}
#' 
#' @examples
#' \dontrun{
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(
#'   arm.data = diabetes.sim, 
#'   varname.t = "Treatment", 
#'   varname.s = "Study"
#' )
#' 
#' # Random effects, consistency model
#' # Binomial family, cloglog link (Hazard Ratio)
#' diabetes.re.c <- nma.model(
#'   data = diabetes.slr,
#'   outcome = "diabetes", 
#'   N = "n",
#'   reference = "Placebo",
#'   family = "binomial",
#'   link = "cloglog",
#'   effects = "random",
#'   type = "consistency",
#'   time = "followup"
#' )
#'  
#' # Run with Stan instead of JAGS
#' diabetes.re.c.stan <- nma.run.stan(
#'   model = diabetes.re.c,
#'   iter_warmup = 1000,
#'   iter_sampling = 2000,
#'   chains = 4
#' )
#' 
#' # Results can be used with all downstream BUGSnet functions
#' nma.league(diabetes.re.c.stan)
#' nma.rank(diabetes.re.c.stan)
#' }
#' 
#' @export
nma.run.stan <- function(
  model,
  monitor = "DEFAULT",
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = chains,
  thin = 1,
  adapt_delta = 0.95,
  max_treedepth = 15,
  seed = NULL,
  refresh = 100,
  show_messages = FALSE,
  ...
) {
  
  # Validate inputs
  if (!inherits(model, 'BUGSnetModel')) {
    stop("'model' must be a valid BUGSnetModel object created using the nma.model function.",
         call. = FALSE)
  }
  
  # Check cmdstanr installation
  check_cmdstan_installation()
  
  # Select appropriate Stan model
  message("Selecting Stan model for: ", model$family, "-", model$link, 
          " with ", model$effects, " effects (", model$type, ")")
  
  stan_file <- select_stan_model(
    family = model$family,
    link = model$link,
    effects = model$effects,
    type = model$type
  )
  
  # Convert data to Stan format
  message("Preparing data for Stan...")
  stan_data <- convert_bugs_to_stan_data(model)
  
  # Compile model
  message("Compiling Stan model...")
  stan_model <- compile_stan_model(stan_file)
  
  if (!show_messages) {
    refresh_compile <- 0
  } else {
    refresh_compile <- refresh
  }
  
  # Run sampling
  message("Starting MCMC sampling with ", chains, " chains...")
  message("  Warmup iterations: ", iter_warmup)
  message("  Sampling iterations: ", iter_sampling)
  message("  Total iterations per chain: ", iter_warmup + iter_sampling)
  
  fit <- tryCatch({
    stan_model$sample(
      data = stan_data,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      parallel_chains = parallel_chains,
      thin = thin,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      seed = seed,
      refresh = refresh,
      show_messages = show_messages,
      ...
    )
  }, error = function(e) {
    stop("Stan sampling failed: ", e$message, call. = FALSE)
  })
  
  message("Sampling completed!")
  
  # Check diagnostics
  message("Checking MCMC diagnostics...")
  sampling_args <- list(
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  diagnostics <- get_stan_diagnostics(fit)
  print_stan_warnings(diagnostics, sampling_args)
  
  # Convert to BUGSnetRun format
  message("Converting output to BUGSnet format...")
  brun <- convert_stan_to_bugsnet(fit, model)
  
  message("Done! Use nma.fit(), nma.league(), nma.rank(), etc. to analyze results.")
  
  return(brun)
}

