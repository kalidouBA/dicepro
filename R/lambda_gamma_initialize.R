# =============================================================================
# lambda_gamma_initialize.R
#
# Public API  : create_gamma_lambda_plot()
#
# NOTE: .custom_space()        is defined in optimisation.R
#       .sample_from_space()   is defined in hypersearch.R
#       globalVariables("bound_label") is declared in utils.R
# =============================================================================


# -----------------------------------------------------------------------------
# create_gamma_lambda_plot  [public]
# -----------------------------------------------------------------------------

#' Visualise the gamma–lambda hyperparameter search space
#'
#' Generates a ggplot2 scatter plot of \eqn{\gamma} vs \eqn{\lambda}
#' sampled from the search space defined by \code{hspaceTechniqueChoose},
#' with feasibility bounds overlaid for \code{"restrictionEspace"}.
#'
#' @param hspaceTechniqueChoose Character scalar.
#'   \describe{
#'     \item{\code{"all"}}{Full independent log-uniform grid for
#'       \code{lambda_}, \code{gamma}, \code{p_prime}. No constraint lines.}
#'     \item{\code{"restrictionEspace"}}{Restricted space with
#'       \code{lambda_ = gamma * lambda_factor},
#'       \code{lambda_factor} \eqn{\in [2, 100]}.
#'       Lower (\eqn{\lambda = 2\gamma}) and upper
#'       (\eqn{\lambda = 100\gamma}) bounds are drawn.}
#'   }
#' @param n_samples Positive integer. Configurations to draw (default \code{200L}).
#'
#' @return A \code{ggplot2} figure object.
#'
#' @export
#' @examples
#' \dontrun{
#' create_gamma_lambda_plot(hspaceTechniqueChoose = "all")
#' create_gamma_lambda_plot(hspaceTechniqueChoose = "restrictionEspace")
#' }
create_gamma_lambda_plot <- function(
    hspaceTechniqueChoose = c("all", "restrictionEspace"),
    n_samples             = 200L) {

  hspaceTechniqueChoose <- match.arg(hspaceTechniqueChoose)

  raw_space <- switch(
    hspaceTechniqueChoose,
    all = list(
      lambda_ = c("loguniform", 1,    1e8),
      gamma   = c("loguniform", 1,    1e8),
      p_prime = c("loguniform", 1e-6, 1)
    ),
    restrictionEspace = .custom_space()   # defined in optimisation.R
  )

  # .sample_from_space() is defined in hypersearch.R
  samples <- lapply(seq_len(n_samples), function(i) .sample_from_space(raw_space))

  plot_df <- data.frame(
    gamma   = vapply(samples, function(s) as.numeric(s$gamma),   numeric(1L)),
    lambda_ = vapply(samples, function(s) as.numeric(s$lambda_), numeric(1L))
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = gamma, y = lambda_)) +
    ggplot2::geom_point(colour = "royalblue", size = 1.5, alpha = 0.45) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  if (hspaceTechniqueChoose == "restrictionEspace") {

    gamma_range <- exp(seq(
      log(min(plot_df$gamma)),
      log(max(plot_df$gamma)),
      length.out = 100L
    ))

    bounds_df <- data.frame(
      gamma       = rep(gamma_range, 2L),
      lambda_     = c(2 * gamma_range, 100 * gamma_range),
      bound_label = rep(
        c("lambda = 2 * gamma (lower bound)",
          "lambda = 100 * gamma (upper bound)"),
        each = 100L
      )
    )

    p <- p +
      ggplot2::geom_line(
        data      = bounds_df,
        mapping   = ggplot2::aes(x = gamma, y = lambda_, colour = bound_label),
        linetype  = "dashed",
        linewidth = 0.8
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          "lambda = 2 * gamma (lower bound)"   = "red",
          "lambda = 100 * gamma (upper bound)" = "green4"
        ),
        name = NULL
      )

    title_str <- paste0(
      "gamma vs lambda \u2014 restrictionEspace ",
      "(lambda = gamma \u00d7 lambda_factor, factor \u2208 [2, 100])"
    )

  } else {
    title_str <- "gamma vs lambda \u2014 all (independent log-uniform sampling)"
  }

  p +
    ggplot2::labs(
      title = title_str,
      x     = expression(gamma ~ "(log scale)"),
      y     = expression(lambda ~ "(log scale)")
    ) +
    ggplot2::theme_bw()
}
