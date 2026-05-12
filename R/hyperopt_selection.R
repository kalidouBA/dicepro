# =============================================================================
# hyperopt_selection.R
#
# Public API  : best_hyperParams()
# Private fns : .select_knee_pareto()
#               .plot_pareto_kraljic()
#               .make_preddata_fname()
#
# Pipeline:
#   1. Filter valid hyper-parameter trials  (abs_constraint <= threshold)
#   2. Compute constraint deviation         (abs_constraint = |1 - constraint|)
#   3. Extract Pareto frontier via KraljicMatrix::get_frontier
#        objectives : minimize frobNorm (x)  AND  abs_constraint (y)
#        → quadrant = "bottom.left"
#   4. Detect the knee point via perpendicular distance to utopia–nadir line
#   5. Generate a ggplot2 visualization
#
# Columns used consistently throughout the pipeline:
#   x-axis  →  frobNorm
#   y-axis  →  abs_constraint
#
# NOTE: utils::globalVariables(c("frobNorm", "abs_constraint")) must be
#       declared in utils.R to silence R CMD CHECK NSE notes.
# =============================================================================


# -----------------------------------------------------------------------------
# .select_knee_pareto  [private]
# -----------------------------------------------------------------------------

#' Select the knee point on a Pareto frontier
#'
#' Identifies the solution with the maximum perpendicular distance to the line
#' connecting the utopia point \eqn{A = (\min f_1, \min f_2)} and the nadir
#' point \eqn{B = (\max f_1, \max f_2)} in normalized objective space.
#'
#' @param df             data.frame of Pareto-optimal solutions.
#' @param loss_col       Character scalar. Column for objective 1 (Frobenius
#'   norm). Default \code{"frobNorm"}.
#' @param constraint_col Character scalar. Column for objective 2 (constraint
#'   deviation). Default \code{"abs_constraint"}.
#' @param eps            Numeric scalar. Guard against zero-range normalization
#'   (default \code{1e-12}).
#' @param return_scores  Logical. Return diagnostic scores when \code{TRUE}
#'   (default); return only \code{best_row} otherwise.
#' @param verbose        Logical. Print a short summary when \code{TRUE}
#'   (default \code{FALSE}).
#'
#' @return Named list (when \code{return_scores = TRUE}):
#' \describe{
#'   \item{best_row}{One-row data.frame of the selected configuration.}
#'   \item{best_index}{Integer row index within \code{df}.}
#'   \item{knee_score}{Perpendicular-distance scores (higher = stronger knee).}
#'   \item{f1_norm}{Normalized \code{loss_col} values.}
#'   \item{f2_norm}{Normalized \code{constraint_col} values.}
#' }
#' When \code{return_scores = FALSE}, returns \code{best_row} directly.
#'
#' @keywords internal
#' @noRd
.select_knee_pareto <- function(df,
                                loss_col       = "frobNorm",
                                constraint_col = "abs_constraint",
                                eps            = 1e-12,
                                return_scores  = TRUE,
                                verbose        = FALSE) {

  stopifnot(is.data.frame(df))
  if (!loss_col %in% names(df))
    stop("Column not found: ", loss_col, call. = FALSE)
  if (!constraint_col %in% names(df))
    stop("Column not found: ", constraint_col, call. = FALSE)

  norm01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < eps) return(rep(0, length(x)))
    (x - rng[1L]) / diff(rng)
  }

  f1_norm <- norm01(df[[loss_col]])
  f2_norm <- norm01(df[[constraint_col]])

  # Utopia A = (min, min),  Nadir B = (max, max) in normalized space
  A <- c(min(f1_norm), min(f2_norm))
  B <- c(max(f1_norm), max(f2_norm))

  # Perpendicular distance from each point to the line AB
  num <- abs(
    (B[2L] - A[2L]) * f1_norm -
      (B[1L] - A[1L]) * f2_norm +
      B[1L] * A[2L] -
      B[2L] * A[1L]
  )
  den        <- sqrt((B[2L] - A[2L])^2 + (B[1L] - A[1L])^2) + eps
  knee_score <- num / den   # higher -> stronger knee

  best_index <- which.max(knee_score)
  best_row   <- df[best_index, , drop = FALSE]

  if (verbose) {
    cat("\n--- Pareto Knee Selection ---\n")
    cat("Selected index :", best_index, "\n")
    cat(loss_col,       "=", df[[loss_col]][best_index],       "\n")
    cat(constraint_col, "=", df[[constraint_col]][best_index], "\n")
    cat("knee_score     =", knee_score[best_index], "\n")
  }

  if (return_scores) {
    list(
      best_row   = best_row,
      best_index = best_index,
      knee_score = knee_score,
      f1_norm    = f1_norm,
      f2_norm    = f2_norm
    )
  } else {
    best_row
  }
}


# -----------------------------------------------------------------------------
# .plot_pareto_kraljic  [private]
# -----------------------------------------------------------------------------

#' Generate a Pareto frontier ggplot2 figure
#'
#' Scatter plot showing:
#' \itemize{
#'   \item All constraint-filtered configurations (grey dots).
#'   \item The Pareto frontier ordered by \code{frobNorm} (red line +
#'         diamond points).
#'   \item The selected knee point (blue star).
#' }
#' Both axes are on a \eqn{\log_{10}} scale.
#'
#' @param data_all             data.frame. All constraint-filtered trials.
#' @param data_frontier        data.frame. Pareto-frontier rows.
#' @param bestPareto           One-row data.frame. Knee-point solution.
#' @param constraint_threshold Numeric scalar. Shown in the plot subtitle.
#'
#' @return A \code{ggplot2} object.
#' @keywords internal
#' @noRd
.plot_pareto_kraljic <- function(data_all,
                                 data_frontier,
                                 bestPareto,
                                 constraint_threshold = 0.1) {

  # Sort frontier by frobNorm -> monotone connecting line, no zigzag
  frontier_ordered <- data_frontier[order(data_frontier$frobNorm), ]

  ggplot2::ggplot() +
    # All valid trials (background)
    ggplot2::geom_point(
      data    = data_all,
      mapping = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour  = "grey60",
      size    = 2L,
      alpha   = 0.6
    ) +
    # Pareto frontier line
    ggplot2::geom_line(
      data      = frontier_ordered,
      mapping   = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour    = "red",
      linewidth = 0.8
    ) +
    # Pareto frontier points (diamonds)
    ggplot2::geom_point(
      data    = frontier_ordered,
      mapping = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour  = "red",
      shape   = 18L,
      size    = 3L
    ) +
    # Knee / best point (blue star)
    ggplot2::geom_point(
      data    = bestPareto,
      mapping = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour  = "blue",
      shape   = 8L,
      size    = 5L,
      stroke  = 1.2
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title    = "Frobenius Norm vs Constraint Deviation \u2014 Pareto Frontier",
      subtitle = sprintf(
        "Constraint threshold: |1 \u2212 constraint| \u2264 %g",
        constraint_threshold
      ),
      x       = "Frobenius Norm (log\u2081\u2080 scale)",
      y       = "|1 \u2212 Constraint| (log\u2081\u2080 scale)",
      caption = "Red diamonds: Pareto frontier  \u2022  Blue star: knee point"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 11L),
      plot.subtitle = ggplot2::element_text(size = 9L,  colour = "grey40"),
      plot.caption  = ggplot2::element_text(size = 8L,  colour = "grey50")
    )
}


# -----------------------------------------------------------------------------
# .make_preddata_fname  [private]
# -----------------------------------------------------------------------------

#' Build standard hyper-parameter result file names
#'
#' @param lambda Numeric vector of \eqn{\lambda} values.
#' @param gamma  Numeric vector of \eqn{\gamma} values.
#' @return Character vector (same length as \code{lambda}).
#' @keywords internal
#' @noRd
.make_preddata_fname <- function(lambda, gamma) {
  paste0(
    "hyperopt_results_lambda_", round(lambda, 2L),
    "_gamma_",                  round(gamma,  2L),
    ".RData"
  )
}


# -----------------------------------------------------------------------------
# best_hyperParams  [public]
# -----------------------------------------------------------------------------

#' Select optimal hyper-parameters using a Pareto frontier
#'
#' Identifies the best \eqn{(\lambda, \gamma)} pair from a set of
#' optimization trials by:
#' \enumerate{
#'   \item Filtering trials where \eqn{|1 - \text{constraint}|} exceeds
#'         \code{constraint_threshold}.
#'   \item Computing \code{abs_constraint = |1 - constraint|}.
#'   \item Extracting the Pareto frontier (minimize \code{frobNorm} AND
#'         \code{abs_constraint}) via \code{KraljicMatrix::get_frontier}
#'         with \code{quadrant = "bottom.left"}.
#'   \item Selecting the knee point via maximum perpendicular distance to
#'         the utopia-nadir line in normalized objective space.
#' }
#'
#' @param trials_df            data.frame of optimization trials as returned
#'   by \code{\link{research_hyperOpt}} (\code{$trials}).  Must contain at
#'   least: \code{lambda_}, \code{gamma}, \code{frobNorm}, \code{constraint}.
#' @param W                    List of \eqn{W} matrices, one per trial.
#' @param H                    List of \eqn{H} matrices, one per trial.
#' @param savePaths            Character scalar. Root directory (kept for
#'   API compatibility; no file I/O is performed here).
#' @param constraint_threshold Numeric scalar. Maximum allowed
#'   \eqn{|1 - \text{constraint}|}; trials above this threshold are discarded
#'   (default \code{0.1}).
#'
#' @return Named list, or \code{invisible(NULL)} when no valid configuration
#'   survives filtering:
#' \describe{
#'   \item{hyperparameters}{List with \code{lambda} and \code{gamma}.}
#'   \item{metrics}{List with \code{frobNorm}, \code{abs_constraint}, and
#'     \code{constraint} for the selected trial.}
#'   \item{trials}{Filtered trials data.frame (constraint-passing rows only).}
#'   \item{W}{The \eqn{W} matrix of the selected trial.}
#'   \item{H}{The \eqn{H} matrix of the selected trial.}
#'   \item{plot}{A \code{ggplot2} Pareto frontier figure.}
#' }
#'
#' @export
best_hyperParams <- function(trials_df,
                             W,
                             H,
                             savePaths,
                             constraint_threshold = 0.1) {

  # ── Guard: empty input ──────────────────────────────────────────────────────
  if (nrow(trials_df) == 0L) {
    warning("trials_df is empty -- no configurations to evaluate.",
            call. = FALSE)
    return(invisible(NULL))
  }

  # ── Required columns ────────────────────────────────────────────────────────
  required_cols <- c("lambda_", "gamma", "frobNorm", "constraint")
  missing_cols  <- setdiff(required_cols, names(trials_df))
  if (length(missing_cols) > 0L) {
    stop(
      "trials_df is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # ── Derived columns (before any row removal) ────────────────────────────────
  trials_df$trial_index    <- seq_len(nrow(trials_df))
  trials_df$log_lambda     <- log10(trials_df$lambda_   + 1)
  trials_df$log_gamma      <- log10(trials_df$gamma     + 1)
  trials_df$log_frob       <- log10(trials_df$frobNorm)
  trials_df$abs_constraint <- abs(1 - trials_df$constraint)

  # ── Constraint filter ───────────────────────────────────────────────────────
  trials_df <- trials_df[
    trials_df$abs_constraint <= constraint_threshold, , drop = FALSE
  ]

  if (nrow(trials_df) == 0L) {
    warning(sprintf(
      "No trials pass the constraint filter (abs_constraint <= %g).",
      constraint_threshold
    ), call. = FALSE)
    return(invisible(NULL))
  }

  # ── Pareto frontier ─────────────────────────────────────────────────────────
  # Minimize BOTH frobNorm (x) and abs_constraint (y)  ->  bottom.left.
  # NSE columns declared via utils::globalVariables() in utils.R.
  frontier_result <- KraljicMatrix::get_frontier(
    data     = trials_df,
    x        = frobNorm,        # objective 1 : minimize Frobenius norm
    y        = abs_constraint,  # objective 2 : minimize constraint deviation
    quadrant = "bottom.left"    # minimize BOTH
  )
  frontier_df <- trials_df[rownames(frontier_result), , drop = FALSE]

  if (nrow(frontier_df) == 0L) {
    warning("Pareto frontier is empty after extraction.", call. = FALSE)
    return(invisible(NULL))
  }

  # ── Knee-point selection ────────────────────────────────────────────────────
  # Perpendicular distance to the utopia-nadir line (scale-free criterion).
  knee_res   <- .select_knee_pareto(
    df             = frontier_df,
    loss_col       = "frobNorm",
    constraint_col = "abs_constraint",
    verbose        = FALSE
  )
  best_row   <- knee_res$best_row
  best_trial <- best_row$trial_index

  # ── Plot ────────────────────────────────────────────────────────────────────
  plot_res <- .plot_pareto_kraljic(
    data_all             = trials_df,
    data_frontier        = frontier_df,
    bestPareto           = best_row,        # same coordinate space as axes
    constraint_threshold = constraint_threshold
  )

  # ── Output ──────────────────────────────────────────────────────────────────
  list(
    hyperparameters = list(
      lambda = best_row$lambda_,
      gamma  = best_row$gamma
    ),
    metrics = list(
      frobNorm       = best_row$frobNorm,
      abs_constraint = best_row$abs_constraint,
      constraint     = best_row$constraint
    ),
    trials = trials_df,
    W      = W[[best_trial]],
    H      = H[[best_trial]],
    plot   = plot_res
  )
}
