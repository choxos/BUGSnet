#' @title Comparison-Adjusted Funnel Plot for Network Meta-Analysis
#' 
#' @description
#' Draws a comparison-adjusted funnel plot to assess funnel plot asymmetry
#' in network meta-analysis, following the method of Chaimani & Salanti (2012).
#' 
#' @param nma A \code{BUGSnetRun} object produced by running \code{nma.run()} or \code{nma.run.stan()}.
#' @param order Mandatory character vector specifying the order of treatments for adjustment
#' @param pooled Character string: "fixed" or "random" (default: based on model)
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis (default: "Standard Error")
#' @param level Confidence level for funnel (default: 0.95)
#' @param col Color(s) for points
#' @param pch Plotting symbol(s)
#' @param legend Logical; add legend (default: TRUE)
#' @param pos.legend Position of legend (default: "topright")
#' @param backtransf Logical; back-transform results
#' @param ... Additional graphical parameters
#' 
#' @return A data frame with comparison-adjusted estimates (invisibly)
#' 
#' @seealso \code{\link{nma.radial}}
#' 
#' @importFrom graphics plot points legend abline
#' @importFrom stats qnorm
#' 
#' @references
#' Chaimani A & Salanti G (2012): Using network meta-analysis to evaluate the 
#' existence of small-study effects in a network of interventions. 
#' Research Synthesis Methods, 3, 161-76
#' 
#' @examples
#' \dontrun{
#' data(diabetes.sim)
#' 
#' diabetes.slr <- data.prep(arm.data = diabetes.sim, 
#'                            varname.t = "Treatment", 
#'                            varname.s = "Study")
#' 
#' diabetes.re <- nma.model(data = diabetes.slr,
#'                          outcome = "diabetes", 
#'                          N = "n",
#'                          reference = "Placebo",
#'                          family = "binomial",
#'                          link = "logit",
#'                          effects = "random")
#' 
#' diabetes.re.res <- nma.run(model = diabetes.re,
#'                            n.adapt = 1000,
#'                            n.burnin = 1000,
#'                            n.iter = 10000)
#' 
#' # Funnel plot (order treatments by recency or some other criterion)
#' nma.funnel(diabetes.re.res, order = c("Placebo", "Lifestyle", "Rosiglitazone", "Metformin"))
#' }
#' 
#' @export
nma.funnel <- function(nma,
                       order,
                       pooled = ifelse(nma$model$effects == "random", "random", "fixed"),
                       xlab = NULL,
                       ylab = "Standard Error",
                       level = 0.95,
                       col = "black",
                       pch = 1,
                       legend = TRUE,
                       pos.legend = "topright",
                       backtransf = TRUE,
                       ...) {
  
  # Validate inputs
  if (!inherits(nma, 'BUGSnetRun'))
    stop("'nma' must be a valid BUGSnetRun object.")
  
  if (missing(order))
    stop("Argument 'order' must be specified. See help page for details.")
  
  # Extract pairwise comparisons from model data
  pairwise_data <- extract_pairwise_data(nma, order)
  
  if (nrow(pairwise_data) == 0)
    stop("No pairwise comparisons found for the specified order.")
  
  # Calculate comparison-adjusted effects
  samples_matrix <- do.call(rbind, nma$samples) %>% data.frame()
  
  # Get pooled estimates for each comparison
  for (i in 1:nrow(pairwise_data)) {
    trt1 <- pairwise_data$treat1[i]
    trt2 <- pairwise_data$treat2[i]
    
    trt1_idx <- which(nma$trt.key == trt1)
    trt2_idx <- which(nma$trt.key == trt2)
    
    # Network estimate
    if (trt1_idx == 1) {
      nma_effect <- samples_matrix[[paste0("d.", trt2_idx, ".")]]
    } else if (trt2_idx == 1) {
      nma_effect <- -samples_matrix[[paste0("d.", trt1_idx, ".")]]
    } else {
      nma_effect <- samples_matrix[[paste0("d.", trt2_idx, ".")]] - samples_matrix[[paste0("d.", trt1_idx, ".")]]
    }
    
    pairwise_data$TE.nma[i] <- mean(nma_effect)
  }
  
  # Calculate comparison-adjusted effect
  pairwise_data$TE.adj <- pairwise_data$TE - pairwise_data$TE.nma
  
  # Back-transform if needed
  if (backtransf && nma$link %in% c("log", "logit", "cloglog")) {
    pairwise_data$TE.adj <- exp(pairwise_data$TE.adj)
    null_val <- 1
  } else {
    null_val <- 0
  }
  
  # Set x-axis label
  if (is.null(xlab)) {
    xlab <- if (backtransf && nma$link %in% c("log", "logit", "cloglog")) {
      paste("Comparison-adjusted", nma$scale)
    } else {
      paste("Comparison-adjusted log", nma$scale)
    }
  }
  
  # Create funnel plot
  plot(pairwise_data$TE.adj, pairwise_data$seTE,
       xlab = xlab, ylab = ylab,
       ylim = rev(range(pairwise_data$seTE)),
       col = col, pch = pch,
       ...)
  
  # Add funnel
  z_val <- qnorm(1 - (1 - level) / 2)
  se_range <- range(pairwise_data$seTE)
  se_seq <- seq(0, max(se_range), length.out = 100)
  
  abline(v = null_val, lty = 2, col = "gray50")
  lines(null_val + z_val * se_seq, se_seq, lty = 3, col = "gray50")
  lines(null_val - z_val * se_seq, se_seq, lty = 3, col = "gray50")
  
  if (legend) {
    comps <- unique(pairwise_data$comparison)
    legend(pos.legend, legend = comps, col = col, pch = pch, cex = 0.8)
  }
  
  invisible(pairwise_data)
}


#' @title Comparison-Adjusted Radial Plot for Network Meta-Analysis
#' 
#' @description
#' Draws a comparison-adjusted radial plot (Galbraith plot) to assess funnel plot
#' asymmetry in network meta-analysis.
#' 
#' @param nma A \code{BUGSnetRun} object
#' @param order Mandatory character vector specifying the order of treatments
#' @param pooled Character string: "fixed" or "random"
#' @param level Confidence level (default: 0.95)
#' @param col Color(s) for points
#' @param pch Plotting symbol(s)
#' @param legend Logical; add legend (default: FALSE)
#' @param pos.legend Position of legend (default: "topright")
#' @param ... Additional graphical parameters
#' 
#' @return A data frame with comparison-adjusted estimates (invisibly)
#' 
#' @seealso \code{\link{nma.funnel}}
#' 
#' @importFrom graphics plot points legend abline
#' @importFrom stats qnorm
#' 
#' @references
#' Chaimani A & Salanti G (2012): Using network meta-analysis to evaluate the 
#' existence of small-study effects in a network of interventions. 
#' Research Synthesis Methods, 3, 161-76
#' 
#' @examples
#' \dontrun{
#' # See examples in nma.funnel
#' }
#' 
#' @export
nma.radial <- function(nma,
                       order,
                       pooled = ifelse(nma$model$effects == "random", "random", "fixed"),
                       level = 0.95,
                       col = "black",
                       pch = 1,
                       legend = FALSE,
                       pos.legend = "topright",
                       ...) {
  
  # Validate inputs
  if (!inherits(nma, 'BUGSnetRun'))
    stop("'nma' must be a valid BUGSnetRun object.")
  
  if (missing(order))
    stop("Argument 'order' must be specified.")
  
  # Extract pairwise comparisons
  pairwise_data <- extract_pairwise_data(nma, order)
  
  if (nrow(pairwise_data) == 0)
    stop("No pairwise comparisons found.")
  
  # Calculate comparison-adjusted effects (same as funnel)
  samples_matrix <- do.call(rbind, nma$samples) %>% data.frame()
  
  for (i in 1:nrow(pairwise_data)) {
    trt1 <- pairwise_data$treat1[i]
    trt2 <- pairwise_data$treat2[i]
    
    trt1_idx <- which(nma$trt.key == trt1)
    trt2_idx <- which(nma$trt.key == trt2)
    
    if (trt1_idx == 1) {
      nma_effect <- samples_matrix[[paste0("d.", trt2_idx, ".")]]
    } else if (trt2_idx == 1) {
      nma_effect <- -samples_matrix[[paste0("d.", trt1_idx, ".")]]
    } else {
      nma_effect <- samples_matrix[[paste0("d.", trt2_idx, ".")]] - samples_matrix[[paste0("d.", trt1_idx, ".")]]
    }
    
    pairwise_data$TE.nma[i] <- mean(nma_effect)
  }
  
  pairwise_data$TE.adj <- pairwise_data$TE - pairwise_data$TE.nma
  
  # Calculate radial plot coordinates
  pairwise_data$x <- 1 / pairwise_data$seTE  # Precision
  pairwise_data$y <- pairwise_data$TE.adj / pairwise_data$seTE  # Standardized effect
  
  # Create radial plot
  z_val <- qnorm(1 - (1 - level) / 2)
  
  plot(pairwise_data$x, pairwise_data$y,
       xlab = "Inverse of Standard Error",
       ylab = "Standardized Comparison-Adjusted Effect",
       col = col, pch = pch,
       ...)
  
  # Add reference lines
  abline(h = 0, lty = 2, col = "gray50")
  abline(h = c(-z_val, z_val), lty = 3, col = "gray50")
  
  if (legend) {
    comps <- unique(pairwise_data$comparison)
    legend(pos.legend, legend = comps, col = col, pch = pch, cex = 0.8)
  }
  
  invisible(pairwise_data)
}


# Helper function to extract pairwise data
extract_pairwise_data <- function(nma, order) {
  # This is a simplified version - would need full implementation
  # to properly extract all pairwise comparisons from the network
  
  data <- nma$model$data
  result_list <- list()
  counter <- 1
  
  # Get study labels
  study_labels <- rownames(data$t_a)
  if (is.null(study_labels)) {
    study_labels <- paste0("Study_", 1:data$ns_a)
  }
  
  for (study in 1:data$ns_a) {
    study_trts_idx <- data$t_a[study, 1:data$na_a[study]]
    study_trts <- nma$trt.key[study_trts_idx]
    
    # For each pair in the study
    if (length(study_trts) >= 2) {
      for (i in 1:(length(study_trts)-1)) {
        for (j in (i+1):length(study_trts)) {
          trt1 <- study_trts[i]
          trt2 <- study_trts[j]
          
          # Calculate effect and SE (simplified - would need proper calculation)
          if (nma$family == "binomial") {
            r1 <- data$r[study, i]
            n1 <- data$n[study, i]
            r2 <- data$r[study, j]
            n2 <- data$n[study, j]
            
            # Add small continuity correction to avoid log(0)
            if (!is.na(r1) && !is.na(r2) && !is.na(n1) && !is.na(n2)) {
              r1_adj <- r1 + 0.5
              r2_adj <- r2 + 0.5
              n1_adj <- n1 + 1
              n2_adj <- n2 + 1
              
              logOR <- log((r2_adj/(n2_adj-r2_adj))/(r1_adj/(n1_adj-r1_adj)))
              seLogOR <- sqrt(1/r1_adj + 1/(n1_adj-r1_adj) + 1/r2_adj + 1/(n2_adj-r2_adj))
              
              result_list[[counter]] <- data.frame(
                studlab = study_labels[study],
                treat1 = trt1,
                treat2 = trt2,
                comparison = paste(trt2, "vs", trt1),
                TE = logOR,
                seTE = seLogOR,
                TE.nma = NA,
                stringsAsFactors = FALSE
              )
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
  
  if (length(result_list) > 0) {
    result <- do.call(rbind, result_list)
  } else {
    result <- data.frame(studlab = character(), treat1 = character(), 
                        treat2 = character(), comparison = character(),
                        TE = numeric(), seTE = numeric(), TE.nma = numeric(),
                        stringsAsFactors = FALSE)
  }
  
  return(result)
}

