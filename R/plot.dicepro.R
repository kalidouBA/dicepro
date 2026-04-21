# =============================================================================
# Plot methods for dicepro objects
#
# Public API : plot.dicepro(), plot_hyperopt(), plot_hyperopt.dicepro()
# Private    : .param_labels, .loss_label, .get_label(), .make_plot()
# =============================================================================


# -----------------------------------------------------------------------------
# plot.dicepro
# -----------------------------------------------------------------------------

#' Plot cell abundance heatmap and error plot
#'
#' Combines a cell-abundance heatmap with a fold-error plot using
#' \pkg{patchwork}. Requires at least two unique iterations in
#' \code{x$Matrix_prediction}.
#'
#' @param x   A \code{dicepro} object as returned by
#'   \code{best_hyperParams}.
#' @param ... Additional arguments (currently unused; reserved for future use).
#'
#' @return A \pkg{patchwork} figure, or \code{invisible(NULL)} with a warning
#'   when only one iteration is present.
#'
#' @method plot dicepro
#' @importFrom patchwork wrap_plots
#' @seealso \code{heatmap_abundances}, \code{metric_plot}
#' @export
plot.dicepro <- function(x, ...) {

  if (length(unique(x$Matrix_prediction$Iterate)) > 1L) {
    patchwork::wrap_plots(
      heatmap_abundances(x$Matrix_prediction),
      metric_plot(x$performs2plot),
      ncol = 1L
    )
  } else {
    warning("Only one unique iteration found -- no plot generated.")
    invisible(NULL)
  }
}


# ---- Parameter label mapping ------------------------------------------------
.param_labels <- list(
  lambda_      = expression(lambda),
  gamma        = expression(gamma),
  gamma_factor = expression(gamma[f]),
  p_prime      = expression(p*minute)
)

.loss_label <- expression(L[obs])

.get_label <- function(p) {
  lab <- .param_labels[[p]]
  if (is.null(lab)) p else lab
}


# -----------------------------------------------------------------------------
# plot_hyperopt  generic
# -----------------------------------------------------------------------------

#' Plot hyperparameter optimisation report
#'
#' Generic function for plotting the hyperparameter search report stored in a
#' \code{dicepro} object. Dispatches to \code{plot_hyperopt.dicepro}.
#'
#' @param x   An object for which a \code{plot_hyperopt} method exists
#'   (currently only \code{dicepro}).
#' @param ... Arguments passed to the method; see
#'   \code{plot_hyperopt.dicepro} for the full list.
#'
#' @return Whatever the dispatched method returns (a \code{patchwork} figure
#'   for \code{dicepro} objects).
#'
#' @seealso \code{plot_hyperopt.dicepro}
#' @export
plot_hyperopt <- function(x, ...) UseMethod("plot_hyperopt")


# -----------------------------------------------------------------------------
# plot_hyperopt.dicepro  method
# -----------------------------------------------------------------------------

#' Plot hyperparameter optimisation report for a dicepro object
#'
#' Builds a scatter-matrix of all evaluated \eqn{(\lambda, \gamma, p')}
#' configurations, colour-coded by loss value, with violin/bar marginals for
#' the top 5\% of trials. Reads all data from \code{x$trials} -- no file I/O.
#'
#' @param x              A \code{dicepro} object. Trials are read from
#'   \code{x$trials}.
#' @param params         Character vector of hyperparameter column names to
#'   display (e.g. \code{c("lambda_", "gamma", "p_prime")}).
#' @param metric         Character scalar. Column used for point size
#'   (default \code{"loss"}).
#' @param loss_metric    Character scalar. Column used as the loss axis
#'   (default \code{"loss"}).
#' @param loss_behaviour Character scalar. Direction of the loss:
#'   \code{"min"} (default) or \code{"max"}.
#' @param not_log        Character vector of parameter names that should
#'   \emph{not} be log-scaled on their axis (default \code{NULL}).
#' @param categorical    Character vector of categorical parameter names;
#'   these are displayed as bar charts in the marginal row
#'   (default \code{NULL}).
#' @param max_deviation  Numeric scalar. Trials with
#'   \eqn{|loss - mean(loss)|} above this value are excluded as outliers
#'   (default \code{NULL}, no exclusion).
#' @param title          Character scalar. Optional title printed above the
#'   combined figure (default \code{NULL}).
#' @param ...            Currently unused. Reserved for future extensions.
#'
#' @return A \code{patchwork} figure which can be printed, saved with
#'   \code{ggplot2::ggsave()}, or embedded in R Markdown documents.
#'
#' @examples
#' \dontrun{
#' out <- dicepro(reference = BlueCode, bulk = CellMixtures,
#'                methodDeconv = "FARDEEP", hp_max_evals = 50L)
#' plot_hyperopt(out, params = c("lambda_", "gamma", "p_prime"))
#' }
#'
#' @rdname plot_hyperopt
#' @import  ggplot2
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @export
plot_hyperopt.dicepro <- function(x,
                                  params,
                                  metric         = "loss",
                                  loss_metric    = "loss",
                                  loss_behaviour = "min",
                                  not_log        = NULL,
                                  categorical    = NULL,
                                  max_deviation  = NULL,
                                  title          = NULL,
                                  ...) {

  trials <- x$trials
  if (is.null(trials) || nrow(trials) == 0L)
    stop("x$trials is NULL or empty.")

  missing_cols <- setdiff(c(params, loss_metric, metric), names(trials))
  if (length(missing_cols) > 0L)
    stop("Columns not found in x$trials: ", paste(missing_cols, collapse = ", "))

  .scale_01 <- function(v) {
    rng <- range(v, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(v)))
    (v - rng[1L]) / diff(rng)
  }

  loss   <- trials[[loss_metric]]
  scores <- trials[[metric]]
  keep   <- !is.na(scores) & !is.na(loss)
  loss   <- loss[keep]
  scores <- scores[keep]
  trials <- trials[keep, , drop = FALSE]

  values <- lapply(stats::setNames(params, params), function(p) trials[[p]])

  if (!is.null(max_deviation)) {
    ok     <- abs(loss - mean(loss, na.rm = TRUE)) < max_deviation
    loss   <- loss[ok]
    scores <- scores[ok]
    values <- lapply(values, `[`, ok)
  }

  categorical <- categorical %||% character(0L)
  for (p in categorical) values[[p]] <- as.character(values[[p]])

  sorted_idx <- do.call(order, c(
    lapply(values[setdiff(params, categorical)], unlist),
    lapply(values[categorical], unlist),
    list(loss, scores)
  ))
  loss   <- loss[sorted_idx]
  scores <- scores[sorted_idx]
  values <- lapply(values, `[`, sorted_idx)

  scores_scaled <- .scale_01(scores)
  lmaxs <- if (loss_behaviour == "min") {
    loss > min(loss, na.rm = TRUE)
  } else {
    loss < max(loss, na.rm = TRUE)
  }
  top_n         <- max(1L, ceiling(length(scores) * 0.05))
  smaxs         <- order(scores, decreasing = TRUE)[seq_len(top_n)]
  cmaxs         <- rep(NA_real_, length(scores))
  cmaxs[smaxs]  <- .scale_01(scores[smaxs])

  df        <- as.data.frame(values)
  df$loss   <- loss
  df$scores <- scores_scaled
  df$smaxs  <- seq_along(scores) %in% smaxs
  df$lmaxs  <- lmaxs
  df$cmaxs  <- cmaxs
  not_log   <- not_log %||% character(0L)

  # ---- Per-panel plot factory -----------------------------------------------
  .make_plot <- function(x_var, y_var, diag = FALSE) {
    if (diag) {
      p <- ggplot(df, aes(x = .data[[x_var]], y = loss)) +
        geom_point(aes(size = scores, colour = lmaxs), alpha = 0.7) +
        geom_point(data = df[df$smaxs, ],
                   aes(size = scores, fill = cmaxs),
                   shape = 21, colour = "black") +
        scale_size_continuous(range = c(1, 10)) +
        scale_colour_manual(values = c("TRUE" = "orange", "FALSE" = "red")) +
        scale_fill_viridis_c() +
        xlab(.get_label(x_var)) +
        ylab(.loss_label) +
        theme_minimal() +
        theme(legend.position = "none")
      if (!x_var %in% not_log) p <- p + scale_x_log10()
    } else {
      p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
        geom_point(aes(size = scores, colour = loss), alpha = 0.7) +
        geom_point(data = df[df$smaxs, ],
                   aes(size = scores, fill = cmaxs),
                   shape = 21, colour = "black") +
        scale_size_continuous(range = c(1, 10)) +
        scale_colour_viridis_c() +
        scale_fill_viridis_c() +
        xlab(.get_label(x_var)) +
        ylab(.get_label(y_var)) +
        theme_minimal() +
        theme(legend.position = "none")
      if (!x_var %in% not_log) p <- p + scale_x_log10()
      if (!y_var %in% not_log) p <- p + scale_y_log10()
    }
    p
  }

  # ---- Scatter matrix -------------------------------------------------------
  scatter_grobs <- vector("list", length(params)^2L)
  k <- 1L
  for (i in seq_along(params))
    for (j in seq_along(params)) {
      scatter_grobs[[k]] <- .make_plot(params[j], params[i],
                                       diag = params[i] == params[j])
      k <- k + 1L
    }

  # ---- Violin / bar row -----------------------------------------------------
  df_best           <- df[df$smaxs, ]
  enough_for_violin <- nrow(df_best) >= 3L

  violin_grobs <- lapply(params, function(p) {
    p_lab <- .get_label(p)

    if (p %in% categorical) {
      ggplot(df_best, aes(x = .data[[p]])) +
        geom_bar(fill = "forestgreen", alpha = 0.3) +
        xlab(p_lab) +
        ylab("Count") +
        theme_minimal()

    } else {
      vp <- ggplot(df_best, aes(x = 1, y = .data[[p]]))

      if (enough_for_violin) {
        vp <- vp +
          geom_violin(fill = "forestgreen", alpha = 0.3) +
          geom_boxplot(width = 0.1, fill = "orange", alpha = 0.7)
      }

      vp <- vp +
        geom_point(aes(colour = cmaxs),
                   position = position_jitter(width = 0.1, seed = 1L)) +
        scale_colour_viridis_c() +
        xlab(p_lab) +
        ylab("") +
        theme_minimal() +
        theme(axis.text.x = element_blank())

      if (!enough_for_violin)
        vp <- vp + labs(subtitle = "(<3 pts - violin omitted)")

      if (!p %in% not_log) vp <- vp + scale_y_log10()
      vp
    }
  })

  # ---- Combine --------------------------------------------------------------
  scatter_panel <- patchwork::wrap_plots(scatter_grobs,
                                         nrow = length(params),
                                         ncol = length(params))

  violin_panel  <- patchwork::wrap_plots(violin_grobs,
                                         nrow = 1L,
                                         ncol = length(params))

  (scatter_panel / violin_panel) +
    patchwork::plot_layout(heights = c(3, 1)) +
    patchwork::plot_annotation(title = title)
}
