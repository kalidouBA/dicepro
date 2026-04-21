# =============================================================================
# Heatmap_error.R
#
# Public  API : heatmap_abundances(), metric_plot()
# Private fns : .theme_dicepro()
# =============================================================================


# -----------------------------------------------------------------------------
# .theme_dicepro  [private]
# -----------------------------------------------------------------------------

#' Shared ggplot2 theme for dicepro plots
#'
#' A clean, publication-ready theme based on \code{theme_bw}, with blank
#' panel backgrounds, visible axis lines, and consistent text sizing. Used
#' internally by all plotting functions to avoid duplication.
#'
#' @param base_size   Numeric. Base font size in points (default \code{9}).
#' @param base_family Character. Font family (default \code{"Helvetica"}).
#'
#' @return A \code{ggplot2} theme object.
#' @keywords internal
#' @noRd
.theme_dicepro <- function(base_size = 9, base_family = "Helvetica") {

  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # Facet strip labels
      strip.text.x     = ggplot2::element_text(size = base_size),
      strip.text.y     = ggplot2::element_text(size = base_size, angle = 90),
      strip.background = ggplot2::element_blank(),
      # Axis text and ticks
      axis.text.x      = ggplot2::element_text(size = base_size),
      axis.text.y      = ggplot2::element_text(size = base_size, hjust = 1),
      axis.ticks       = ggplot2::element_line(colour = "black"),
      axis.title.x     = ggplot2::element_text(size = base_size),
      axis.title.y     = ggplot2::element_text(size = base_size, angle = 90),
      axis.line.x      = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y      = ggplot2::element_line(color = "black", linewidth = 0.5),
      # Panel and plot backgrounds
      panel.background = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      plot.margin      = ggplot2::unit(c(0.5, 1, 1, 1), "lines"),
      plot.title       = ggplot2::element_text(hjust = 0.5, size = base_size + 1L,
                                               face = "bold"),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = base_size,
                                               colour = "grey40"),
      # Legend
      legend.text      = ggplot2::element_text(size = base_size),
      legend.title     = ggplot2::element_text(size = base_size, face = "bold"),
      legend.key.size  = ggplot2::unit(0.5, "cm"),
      legend.background = ggplot2::element_blank()
    )
}


# -----------------------------------------------------------------------------
# heatmap_abundances  [public]
# -----------------------------------------------------------------------------

#' Heatmap of cell-type abundances
#'
#' Generates a raster heatmap showing the estimated abundance of each cell
#' type across NMF iterations (or samples). Colour encodes abundance using
#' the viridis scale (reversed: dark = high abundance).
#'
#' @param res2plot  A \code{data.frame} with:
#'   \itemize{
#'     \item one column named \code{Iterate} (iteration or sample index,
#'       numeric or integer);
#'     \item one column per cell type containing non-negative numeric
#'       abundance values.
#'   }
#' @param base_size  Numeric. Base font size in points (default \code{9}).
#' @param title      Character. Plot title (default \code{NULL}, no title).
#' @param midpoint   Numeric or \code{NULL}. When not \code{NULL}, a
#'   diverging viridis scale is used with this value as the midpoint,
#'   useful for highlighting deviations from a reference abundance.
#'   Default \code{NULL} (sequential scale).
#'
#' @return A \code{ggplot} object. Print or save with
#'   \code{ggplot2::ggsave()}.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Iterate   = 1:10,
#'   CellTypeA = runif(10),
#'   CellTypeB = runif(10)
#' )
#' heatmap_abundances(df)
#' heatmap_abundances(df, title = "Estimated abundances", midpoint = 0.5)
#' }
#'
#' @export
heatmap_abundances <- function(res2plot,
                               base_size = 9,
                               title     = NULL,
                               midpoint  = NULL) {

  # Silence R CMD CHECK NSE notes (declared in utils.R)
  Iterate <- Cell_Type <- Abundances <- NULL

  # ---- Input validation ----------------------------------------------------
  if (!is.data.frame(res2plot))
    stop("'res2plot' must be a data.frame.", call. = FALSE)
  if (!"Iterate" %in% names(res2plot))
    stop("'res2plot' must contain an 'Iterate' column.", call. = FALSE)

  iter_col <- res2plot$Iterate
  ct_cols  <- setdiff(names(res2plot), "Iterate")

  if (length(ct_cols) == 0L)
    stop("'res2plot' must contain at least one cell-type column.", call. = FALSE)

  # ---- Reshape to long format ----------------------------------------------
  data2plot <- data.frame(
    Iterate    = rep(iter_col, times = length(ct_cols)),
    Cell_Type  = rep(ct_cols,  each  = length(iter_col)),
    Abundances = unlist(res2plot[ct_cols], use.names = FALSE),
    stringsAsFactors = FALSE
  )

  # Preserve cell-type column order from the original data.frame
  data2plot$Cell_Type <- factor(data2plot$Cell_Type, levels = rev(ct_cols))

  # ---- Build plot ----------------------------------------------------------
  p <- ggplot2::ggplot(
    data2plot,
    ggplot2::aes(x = Iterate, y = Cell_Type, fill = Abundances)
  ) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::labs(
      x     = "Iteration",
      y     = "Cell type",
      fill  = "Abundance",
      title = title
    ) +
    .theme_dicepro(base_size = base_size) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1,
                                           size = base_size),
      legend.position = "right"
    )

  # ---- Colour scale --------------------------------------------------------
  if (is.null(midpoint)) {
    p <- p + ggplot2::scale_fill_viridis_c(direction = -1, option = "viridis")
  } else {
    # Diverging scale centred on midpoint
    p <- p + ggplot2::scale_fill_gradient2(
      low      = "#440154",
      mid      = "#21908C",
      high     = "#FDE725",
      midpoint = midpoint
    )
  }

  p
}


# -----------------------------------------------------------------------------
# metric_plot  [public]
# -----------------------------------------------------------------------------

#' Performance metric line plot across iterations
#'
#' Plots the evolution of one or more scalar performance metrics (e.g.,
#' NRMSE, R\eqn{^2}, loss) over NMF iterations. Useful for monitoring
#' convergence and diagnosing early stopping.
#'
#' @param perf2plot  A \code{data.frame} with at least two columns:
#'   \itemize{
#'     \item \code{Iterate} — iteration index (numeric or integer).
#'     \item \code{metric}  — scalar performance value at each iteration.
#'   }
#'   Optionally a \code{group} column (character or factor) to draw
#'   multiple coloured lines, one per group.
#' @param ylab       Character. Y-axis label
#'   (default \code{"Error between folds"}).
#' @param title      Character. Plot title (default \code{NULL}).
#' @param base_size  Numeric. Base font size in points (default \code{9}).
#' @param add_smooth Logical. When \code{TRUE}, overlays a LOESS smoothing
#'   line. Useful for noisy convergence curves (default \code{FALSE}).
#' @param show_points Logical. When \code{TRUE}, adds individual data points
#'   on top of the line (default \code{FALSE}).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Iterate = 1:50,
#'   metric  = cumsum(runif(50, -0.01, 0.1))
#' )
#' metric_plot(df)
#' metric_plot(df, ylab = "NRMSE", add_smooth = TRUE)
#'
#' # Multi-group
#' df2 <- data.frame(
#'   Iterate = rep(1:20, 2),
#'   metric  = c(cumsum(runif(20, 0, 0.1)), cumsum(runif(20, 0, 0.05))),
#'   group   = rep(c("Fold 1", "Fold 2"), each = 20)
#' )
#' metric_plot(df2)
#' }
#'
#' @export
metric_plot <- function(perf2plot,
                        ylab         = "Error between folds",
                        title        = NULL,
                        base_size    = 9,
                        add_smooth   = FALSE,
                        show_points  = FALSE) {

  # Silence R CMD CHECK NSE notes
  Iterate <- metric <- group <- NULL

  # ---- Input validation ----------------------------------------------------
  if (!is.data.frame(perf2plot))
    stop("'perf2plot' must be a data.frame.", call. = FALSE)

  required_cols <- c("Iterate", "metric")
  missing_cols  <- setdiff(required_cols, names(perf2plot))
  if (length(missing_cols) > 0L)
    stop(sprintf(
      "'perf2plot' must contain columns: %s. Missing: %s.",
      paste(required_cols, collapse = ", "),
      paste(missing_cols,  collapse = ", ")
    ), call. = FALSE)

  has_group <- "group" %in% names(perf2plot)

  # ---- Base mapping --------------------------------------------------------
  base_aes <- if (has_group) {
    ggplot2::aes(x = Iterate, y = metric, colour = group, group = group)
  } else {
    ggplot2::aes(x = Iterate, y = metric)
  }

  p <- ggplot2::ggplot(perf2plot, base_aes)

  # ---- Layers --------------------------------------------------------------
  if (has_group) {
    p <- p + ggplot2::geom_line(linewidth = 0.7, alpha = 0.9)
  } else {
    p <- p + ggplot2::geom_line(colour = "black", linewidth = 0.7)
  }

  if (show_points) {
    if (has_group) {
      p <- p + ggplot2::geom_point(size = 1.5, alpha = 0.8)
    } else {
      p <- p + ggplot2::geom_point(colour = "black", size = 1.5)
    }
  }

  if (add_smooth) {
    smooth_aes <- if (has_group) {
      ggplot2::aes(x = Iterate, y = metric, colour = group)
    } else {
      ggplot2::aes(x = Iterate, y = metric)
    }
    p <- p + ggplot2::geom_smooth(
      mapping = smooth_aes,
      method  = "loess",
      formula = y ~ x,
      se      = FALSE,
      linewidth = 1.2,
      linetype  = "dashed",
      alpha     = 0.6
    )
  }

  # ---- Labels and theme ----------------------------------------------------
  p <- p +
    ggplot2::labs(
      x      = "Iteration",
      y      = ylab,
      title  = title,
      colour = if (has_group) "Group" else NULL
    ) +
    .theme_dicepro(base_size = base_size)

  if (!has_group) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = "right")
  }

  p
}
