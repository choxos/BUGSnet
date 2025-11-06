#' @title Forest Plot Generic
#' 
#' @description
#' Generic function for forest plots
#' 
#' @param x An object for which a forest plot method exists
#' @param ... Additional arguments passed to methods
#' 
#' @export
forest <- function(x, ...) {
  UseMethod("forest")
}


#' @title Forest Plot for Direct and Indirect Evidence
#' 
#' @description
#' Creates a forest plot showing network meta-analysis estimates alongside
#' direct and indirect evidence estimates, styled similar to netmeta package.
#' 
#' @param x An object of class \code{nma.netsplit} produced by \code{nma.netsplit()}.
#' @param show Character string indicating which comparisons to show: "all", "with.direct" (default), 
#'   "both" (comparisons with both direct and indirect evidence)
#' @param overall Logical indicating whether to show network estimates (default: TRUE)
#' @param direct Logical indicating whether to show direct estimates (default: TRUE)
#' @param indirect Logical indicating whether to show indirect estimates (default: TRUE)
#' @param only.reference Logical indicating whether to show only comparisons with reference group
#' @param sortvar Optional vector to sort comparisons
#' @param backtransf Logical indicating whether to back-transform results (default: TRUE for OR, RR, HR)
#' @param digits Number of digits for estimates (default: 2)
#' @param col.diamond Color for diamond (network estimate)
#' @param col.square Color for squares (direct/indirect estimates)
#' @param equal.size Logical; if TRUE, all estimates have same size (default: TRUE)
#' @param text.overall Label for network estimates (default: "Network")
#' @param text.direct Label for direct estimates (default: "Direct")
#' @param text.indirect Label for indirect estimates (default: "Indirect")
#' @param xlab X-axis label
#' @param ... Additional graphical parameters
#' 
#' @return A ggplot2 object with the forest plot
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_vline facet_grid theme_bw 
#'   theme element_text element_blank labs scale_shape_manual scale_color_manual coord_cartesian
#' @importFrom dplyr %>% filter mutate arrange
#' 
#' @seealso \code{\link{nma.netsplit}}
#' 
#' @examples
#' \dontrun{
#' # See examples in nma.netsplit
#' }
#' 
#' @method forest nma.netsplit
#' @export
forest.nma.netsplit <- function(x,
                                show = "with.direct",
                                overall = TRUE,
                                direct = TRUE,
                                indirect = TRUE,
                                only.reference = FALSE,
                                sortvar = NULL,
                                backtransf = TRUE,
                                digits = 2,
                                col.diamond = "gray30",
                                col.square = "gray60",
                                equal.size = TRUE,
                                text.overall = "Network",
                                text.direct = "Direct",
                                text.indirect = "Indirect",
                                xlab = NULL,
                                ...) {
  
  # Validate show argument
  show <- match.arg(show, c("all", "with.direct", "both"))
  
  # Determine which comparisons to include
  if (show == "all") {
    sel <- rep(TRUE, length(x$comparison))
  } else if (show == "with.direct") {
    sel <- x$k > 0
  } else if (show == "both") {
    sel <- x$k > 0 & !is.na(x$indirect$TE)
  }
  
  if (only.reference && !is.null(x$reference.group)) {
    ref_sel <- grepl(x$reference.group, x$comparison, fixed = TRUE)
    sel <- sel & ref_sel
  }
  
  # Build data frame for plotting
  plot_data <- data.frame()
  
  for (i in which(sel)) {
    comp <- x$comparison[i]
    
    if (overall) {
      plot_data <- rbind(plot_data, data.frame(
        comparison = comp,
        type = text.overall,
        TE = x$nma$TE[i],
        lower = x$nma$lower[i],
        upper = x$nma$upper[i],
        se = x$nma$seTE[i],
        shape = "diamond",
        color = col.diamond,
        stringsAsFactors = FALSE
      ))
    }
    
    if (direct && x$k[i] > 0) {
      plot_data <- rbind(plot_data, data.frame(
        comparison = comp,
        type = text.direct,
        TE = x$direct$TE[i],
        lower = x$direct$lower[i],
        upper = x$direct$upper[i],
        se = x$direct$seTE[i],
        shape = "square",
        color = col.square,
        stringsAsFactors = FALSE
      ))
    }
    
    if (indirect && !is.na(x$indirect$TE[i])) {
      plot_data <- rbind(plot_data, data.frame(
        comparison = comp,
        type = text.indirect,
        TE = x$indirect$TE[i],
        lower = x$indirect$lower[i],
        upper = x$indirect$upper[i],
        se = x$indirect$seTE[i],
        shape = "square",
        color = col.square,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(plot_data) == 0) {
    stop("No data to plot with the specified options")
  }
  
  # Back-transform if needed
  if (backtransf && x$link %in% c("log", "logit", "cloglog")) {
    plot_data$TE <- exp(plot_data$TE)
    plot_data$lower <- exp(plot_data$lower)
    plot_data$upper <- exp(plot_data$upper)
    null_value <- 1
  } else {
    null_value <- 0
  }
  
  # Set x-axis label
  if (is.null(xlab)) {
    xlab <- if (backtransf && x$link %in% c("log", "logit", "cloglog")) {
      x$scale
    } else {
      paste("log", x$scale)
    }
  }
  
  # Create forest plot
  p <- ggplot(plot_data, aes(x = TE, y = comparison)) +
    geom_vline(xintercept = null_value, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), 
                   height = 0.3, size = 0.5) +
    geom_point(aes(shape = shape, fill = color), 
               size = if(equal.size) 3 else (1/plot_data$se * 5)) +
    facet_grid(type ~ ., scales = "free_y", space = "free_y") +
    labs(x = xlab, y = "") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # Add scale transformations
  if (backtransf && x$link %in% c("log", "logit", "cloglog")) {
    p <- p + scale_x_log10()
  }
  
  return(p)
}

