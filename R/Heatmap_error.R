# =============================================================================
# Visualization Functions — Cell Abundances & Performance Metrics
#
# Public  API : heatmap_abundances(), metric_plot()
# Private fns : .theme_dicepro()
# -----------------------------------------------------------------------------
# .theme_dicepro  [private]
# -----------------------------------------------------------------------------

#' Shared ggplot2 Theme for dicepro Plots
#'
#' A clean, publication-ready theme based on \code{theme_bw}, with blank panels,
#' visible axis lines, and consistent text sizing. Used internally by all
#' plotting functions to avoid duplication.
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
      strip.text.x     = ggplot2::element_text(size = base_size),
      strip.text.y     = ggplot2::element_text(size = base_size, angle = 90),
      strip.background = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(size = base_size),
      axis.text.y      = ggplot2::element_text(size = base_size, hjust = 1),
      axis.ticks       = ggplot2::element_line(colour = "black"),
      axis.title.x     = ggplot2::element_text(size = base_size),
      axis.title.y     = ggplot2::element_text(size = base_size, angle = 90),
      axis.line.x      = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y      = ggplot2::element_line(color = "black", linewidth = 0.5),
      panel.background = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background  = ggplot2::element_blank(),
      plot.margin      = ggplot2::unit(c(0.5, 1, 1, 1), "lines"),
      plot.title       = ggplot2::element_text(hjust = 0.5),
      legend.text      = ggplot2::element_text(size = base_size),
      legend.key.size  = ggplot2::unit(0.5, "cm")
    )
}


# -----------------------------------------------------------------------------
# heatmap_abundances  [public]
# -----------------------------------------------------------------------------

#' Heatmap of Cell Abundances Across Iterations
#'
#' Generates a raster heatmap showing the abundance of each cell type at every
#' NMF iteration. Colour encodes abundance level using the viridis scale
#' (reversed: dark = high abundance).
#'
#' @param res2plot  A data.frame with one column named \code{Iterate} (iteration
#'   index) and one column per cell type containing numeric abundance values.
#' @param base_size Numeric. Base font size in points (default \code{9}).
#'
#' @return A \code{ggplot} object. Print it or save it with
#'   \code{ggplot2::ggsave()}.
#'
#' @examples
#' \dontrun{
# df <- data.frame(
#   Iterate  = 1:10,
#   CellTypeA = runif(10),
#   CellTypeB = runif(10)
# )
#' heatmap_abundances(df)
#' }
#'
#' @export
heatmap_abundances <- function(res2plot, base_size = 9) {

  Iterate <- Cell_Type <- Abundances <- NULL
  iter_col   <- res2plot$Iterate
  ct_cols    <- setdiff(names(res2plot), "Iterate")

  data2plot <- data.frame(
    Iterate    = rep(iter_col, times = length(ct_cols)),
    Cell_Type  = rep(ct_cols,  each  = length(iter_col)),
    Abundances = unlist(res2plot[ct_cols], use.names = FALSE),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    data2plot,
    ggplot2::aes(x = Iterate, y = Cell_Type, fill = Abundances)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(direction = -1) +
    ggplot2::labs(x = "Iteration", y = "Cell Type", fill = "Abundance") +
    .theme_dicepro(base_size = base_size)
}


# -----------------------------------------------------------------------------
# metric_plot  [public]
# -----------------------------------------------------------------------------

#' Performance Metric Line Plot Across Iterations
#'
#' Plots the evolution of a scalar performance metric (e.g., NRMSE, R²) over
#' NMF iterations. Useful for monitoring convergence and diagnosing early
#' stopping.
#'
#' @param perf2plot A data.frame with at least two columns:
#'   \itemize{
#'     \item \code{Iterate} — iteration index (numeric or integer).
#'     \item \code{metric}  — scalar performance value at each iteration.
#'   }
#' @param ylab      Character. Y-axis label (default \code{"Error between folds"}).
#' @param base_size Numeric. Base font size in points (default \code{9}).
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
#' metric_plot(df, ylab = "NRMSE")
#' }
#'
#' @export
metric_plot <- function(perf2plot,
                        ylab      = "Error between folds",
                        base_size = 9) {
  Iterate <- metric <- NULL

  ggplot2::ggplot(
    perf2plot,
    ggplot2::aes(x = Iterate, y = metric)
  ) +
    ggplot2::geom_line(colour = "black") +
    ggplot2::labs(x = "Iteration", y = ylab) +
    .theme_dicepro(base_size = base_size)
}
