#' @title Heat Plot for Network Meta-Analysis
#' 
#' @description
#' Produces a heat plot containing treatment estimates with confidence intervals
#' for all possible pairwise comparisons, where a color scale represents the values
#' of relative treatment effects. Similar to a league table but with color coding.
#' 
#' @param nma A \code{BUGSnetRun} object produced by running \code{nma.run()} or \code{nma.run.stan()}.
#' @param central.tdcy The posterior statistic: "mean" or "median" (default: "median")
#' @param order Optional vector of treatment names specifying the order
#' @param low.colour Color for low treatment effects (default: "red")
#' @param mid.colour Color for null treatment effects (default: "white")
#' @param high.colour Color for high treatment effects (default: "springgreen4")
#' @param size Text size for estimates (default: 4)
#' @param size.trt Text size for treatment names (default: 12)
#' @param digits Number of digits for estimates (default: 2)
#' @param backtransf Logical; back-transform results (default: TRUE for OR, RR, HR)
#' @param cov.value For meta-regression, the covariate value (default: NULL)
#' 
#' @return A ggplot2 object with the heat plot
#' 
#' @seealso \code{\link{nma.league}}
#' 
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 theme_minimal 
#'   theme element_text element_blank scale_x_discrete scale_y_discrete labs
#' @importFrom dplyr %>% summarise_all
#' @importFrom tidyr gather
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
#' # Create heat plot
#' nma.heatplot(diabetes.re.res)
#' }
#' 
#' @export
nma.heatplot <- function(nma,
                        central.tdcy = "median",
                        order = NULL,
                        low.colour = "red",
                        mid.colour = "white",
                        high.colour = "springgreen4",
                        size = 4,
                        size.trt = 12,
                        digits = 2,
                        backtransf = TRUE,
                        cov.value = NULL) {
  
  # Validate inputs
  if (!inherits(nma, 'BUGSnetRun'))
    stop("'nma' must be a valid BUGSnetRun object created using nma.run() or nma.run.stan().")
  
  if (nma$model$type == "inconsistency")
    stop("heatplot is not yet implemented for inconsistency models.")
  
  if (!is.null(nma$model$covariate) && is.null(cov.value))
    stop("cov.value must be specified for meta-regression")
  
  # Extract treatment information
  trts <- nma$trt.key
  n.trts <- length(trts)
  
  # Extract samples
  samples_matrix <- do.call(rbind, nma$samples) %>% data.frame()
  d_samples <- samples_matrix %>% select(starts_with("d."))
  colnames(d_samples) <- trts
  
  # Add meta-regression if applicable
  if (!is.null(cov.value)) {
    beta_samples <- samples_matrix %>% select(starts_with("beta."))
    colnames(beta_samples) <- trts
    d_samples <- d_samples + beta_samples * (cov.value - nma$model$mean.cov)
  }
  
  # Calculate pairwise comparisons
  comp_matrix <- matrix(NA, nrow = n.trts, ncol = n.trts)
  lower_matrix <- matrix(NA, nrow = n.trts, ncol = n.trts)
  upper_matrix <- matrix(NA, nrow = n.trts, ncol = n.trts)
  
  for (i in 1:n.trts) {
    for (j in 1:n.trts) {
      if (i != j) {
        # Calculate relative effect (row vs column)
        rel_effect <- d_samples[, j] - d_samples[, i]
        
        if (central.tdcy == "mean") {
          comp_matrix[i, j] <- mean(rel_effect)
        } else if (central.tdcy == "median") {
          comp_matrix[i, j] <- median(rel_effect)
        } else {
          stop("central.tdcy must be 'mean' or 'median'")
        }
        
        lower_matrix[i, j] <- quantile(rel_effect, 0.025)
        upper_matrix[i, j] <- quantile(rel_effect, 0.975)
      }
    }
  }
  
  rownames(comp_matrix) <- colnames(comp_matrix) <- trts
  rownames(lower_matrix) <- colnames(lower_matrix) <- trts
  rownames(upper_matrix) <- colnames(upper_matrix) <- trts
  
  # Back-transform if needed
  if (backtransf && nma$link != "identity") {
    comp_matrix <- exp(comp_matrix)
    lower_matrix <- exp(lower_matrix)
    upper_matrix <- exp(upper_matrix)
    null_value <- 1
  } else {
    null_value <- 0
  }
  
  # Reorder if specified
  if (!is.null(order)) {
    order_idx <- match(order, trts)
    order_idx <- order_idx[!is.na(order_idx)]
    if (length(order_idx) > 0) {
      comp_matrix <- comp_matrix[order_idx, order_idx]
      lower_matrix <- lower_matrix[order_idx, order_idx]
      upper_matrix <- upper_matrix[order_idx, order_idx]
      trts <- trts[order_idx]
    }
  }
  
  # Prepare data for ggplot
  plot_data <- data.frame()
  for (i in 1:nrow(comp_matrix)) {
    for (j in 1:ncol(comp_matrix)) {
      if (i != j) {
        plot_data <- rbind(plot_data, data.frame(
          Treatment = trts[i],
          Comparator = trts[j],
          TE = comp_matrix[i, j],
          lower = lower_matrix[i, j],
          upper = upper_matrix[i, j],
          label = sprintf("%.*f\n[%.*f, %.*f]", 
                         digits, comp_matrix[i, j],
                         digits, lower_matrix[i, j],
                         digits, upper_matrix[i, j]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Create heat plot
  Treatment <- Comparator <- TE <- label <- NULL  # Avoid R CMD check notes
  
  p <- ggplot(plot_data, aes(x = Comparator, y = Treatment, fill = TE)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = label), size = size, color = "black") +
    scale_fill_gradient2(
      low = low.colour,
      mid = mid.colour,
      high = high.colour,
      midpoint = null_value,
      name = nma$scale
    ) +
    scale_x_discrete(position = "top", limits = trts) +
    scale_y_discrete(limits = rev(trts)) +
    labs(
      title = paste("Network Meta-Analysis Heat Plot"),
      subtitle = paste(nma$scale, "-", nma$model$effects, "effects model"),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, size = size.trt, face = "bold"),
      axis.text.y = element_text(size = size.trt, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  
  return(p)
}

