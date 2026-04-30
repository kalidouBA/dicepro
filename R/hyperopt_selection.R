# =============================================================================
# hyperopt_selection.R
#
# Public API  : best_hyperParams()
# Private fns : .select_knee_pareto()
#               .plot_pareto_kraljic()
#               .make_preddata_fname()
#
# Pipeline:
#   1. Filter valid hyper-parameter trials
#   2. Compute constraint deviation
#   3. Extract Pareto frontier (KraljicMatrix::get_frontier)
#   4. Detect the knee point (best compromise)
#   5. Generate a ggplot2 visualization
#
# NOTE: globalVariables(c("frob_H", "loss")) declared in utils.R
# =============================================================================


# -----------------------------------------------------------------------------
# .select_knee_pareto  [private]
# -----------------------------------------------------------------------------

#' Select the knee point on a Pareto frontier
#'
#' Identifies the configuration that maximises the perpendicular distance to
#' the line connecting the utopia \eqn{(0,0)} and nadir \eqn{(1,1)} points
#' in a bi-objective space. Both objectives are min-max normalised before
#' scoring, making the procedure fully scale-free.
#'
#' @param df             data.frame of Pareto-optimal solutions.
#' @param loss_col       Character scalar. Column name of the loss objective.
#' @param constraint_col Character scalar. Column name of the constraint
#'   deviation objective.
#' @param eps            Numeric scalar. Guard against zero-range objectives
#'   (default \code{1e-12}).
#' @param return_scores  Logical. Return diagnostic scores when \code{TRUE}
#'   (default); return only \code{best_row} otherwise.
#'
#' @return Named list (when \code{return_scores = TRUE}):
#' \describe{
#'   \item{best_row}{One-row data.frame of the selected configuration.}
#'   \item{best_index}{Integer row index within \code{df}.}
#'   \item{knee_score}{Numeric vector of perpendicular-distance scores.}
#'   \item{f1_norm}{Normalised loss values.}
#'   \item{f2_norm}{Normalised constraint-deviation values.}
#' }
#' When \code{return_scores = FALSE}, returns \code{best_row} directly.
#'
#' @keywords internal
#' @noRd
.select_knee_pareto <- function(df,
                                loss_col       = "loss",
                                constraint_col = "abs_constraint",
                                eps            = 1e-12,
                                return_scores  = TRUE) {

  stopifnot(is.data.frame(df))
  if (!loss_col       %in% names(df))
    stop("Column not found: ", loss_col,       call. = FALSE)
  if (!constraint_col %in% names(df))
    stop("Column not found: ", constraint_col, call. = FALSE)

  norm01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < eps) return(rep(0, length(x)))
    (x - rng[1L]) / diff(rng)
  }

  f1_norm    <- norm01(df[[loss_col]])
  f2_norm    <- norm01(df[[constraint_col]])
  knee_score <- abs(f1_norm - f2_norm) / sqrt(2)
  best_index <- which.max(knee_score)
  best_row   <- df[best_index, , drop = FALSE]

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
#' Scatter plot showing all valid configurations (grey), the Pareto frontier
#' (red line + diamonds), and the selected knee point (blue star).
#'
#' @param data_all             data.frame of all valid trials.
#' @param data_frontier        data.frame of Pareto-frontier points.
#' @param bestPareto           One-row data.frame of the selected solution.
#' @param constraint_threshold Numeric scalar. Shown in the plot title.
#'
#' @return A \code{ggplot2} object.
#' @keywords internal
#' @noRd
.plot_pareto_kraljic <- function(data_all,
                                 data_frontier,
                                 bestPareto,
                                 constraint_threshold = 0.1) {

  frontier_ordered <- data_frontier[order(data_frontier$frobNorm), ]

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data    = data_all,
      mapping = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour  = "grey60",
      size    = 2L,
      alpha   = 0.6
    ) +
    ggplot2::geom_line(
      data      = frontier_ordered,
      mapping   = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour    = "red",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data    = frontier_ordered,
      mapping = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour  = "red",
      shape   = 18L,
      size    = 3L
    ) +
    ggplot2::geom_point(
      data    = bestPareto,
      mapping = ggplot2::aes(x = frobNorm, y = abs_constraint),
      colour  = "blue",
      shape   = 8L,
      size    = 5L,
      stroke  = 1.2
    ) +
    ggplot2::scale_x_log10(trans = "reverse") +
    ggplot2::scale_y_log10(trans = "reverse") +
    ggplot2::labs(
      title = "Frobenius Norm vs Constraint Deviation \u2014 Pareto Frontier",
      x     = "Frobenius Norm (log scale reversed)",
      y     = "|1 - Constraint| (log scale reversed)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11L)
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
# best_hyper-Params  [public]
# -----------------------------------------------------------------------------

#' Select optimal hyper-parameters using a Pareto frontier
#'
#' Identifies the best \eqn{(\lambda, \gamma)} pair from a set of
#' optimization trials by:
#' \enumerate{
#'   \item Filtering trials that exceed the constraint threshold.
#'   \item Computing constraint deviation (\code{abs_constraint}).
#'   \item Extracting the Pareto frontier (Frobenius norm vs. loss) via
#'     \code{KraljicMatrix::get_frontier}.
#'   \item Selecting the knee point — the best loss/constraint trade-off.
#' }
#'
#' @param trials_df            data.frame of optimization trials as returned
#'   by \code{\link{research_hyperOpt}} (\code{$trials}).
#' @param W                    List of \eqn{W} matrices, one per trial.
#' @param H                    List of \eqn{H} matrices, one per trial.
#' @param savePaths            Character scalar. Root directory (kept for
#'   API compatibility; no file I/O performed here).
#' @param constraint_threshold Numeric scalar. Maximum allowed
#'   \eqn{|1 - \text{constraint}|}; trials above this are discarded
#'   (default \code{0.1}).
#'
#' @return Named list, or \code{invisible(NULL)} when no valid configuration
#'   survives filtering:
#' \describe{
#'   \item{hyperparameters}{List with \code{lambda} and \code{gamma}.}
#'   \item{metrics}{List with \code{loss} and \code{constraint}.}
#'   \item{trials}{Filtered trials data.frame.}
#'   \item{W}{The \eqn{W} matrix of the selected trial.}
#'   \item{H}{The \eqn{H} matrix of the selected trial (first column
#'     dropped — the \code{Mixture} label column).}
#'   \item{plot}{ggplot2 Pareto frontier figure.}
#' }
#'
#' @export
best_hyperParams <- function(trials_df,
                             W,
                             H,
                             savePaths,
                             constraint_threshold = 0.1) {

  # ---- Guard: empty input --------------------------------------------------
  if (nrow(trials_df) == 0L) {
    warning("trials_df is empty -- no configurations to evaluate.",
            call. = FALSE)
    return(invisible(NULL))
  }

  # ---- Derived columns (attach index before any row removal) ---------------
  trials_df$trial_index    <- seq_len(nrow(trials_df))
  trials_df$log_lambda     <- log10(trials_df$lambda_ + 1)
  trials_df$log_gamma      <- log10(trials_df$gamma   + 1)
  trials_df$log_frob       <- log10(trials_df$frobNorm)
  trials_df$abs_constraint <- abs(1 - trials_df$constraint)

  # ---- Constraint filter ---------------------------------------------------
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

  # ---- Pareto frontier -----------------------------------------------------
  # frob_H and loss are column names of trials_df; declared in utils.R via
  # utils::globalVariables() to silence R CMD CHECK NSE notes.
  frontier_result <- KraljicMatrix::get_frontier(
    data     = trials_df,
    x        = frob_H,
    y        = loss,
    quadrant = "bottom.right"
  )
  frontier_df <- trials_df[rownames(frontier_result), , drop = FALSE]

  # ---- Knee-point selection ------------------------------------------------
  # Uses frob_H (x-axis) and var_H (y-axis) as the two objectives.
  knee_res   <- .select_knee_pareto(
    df             = frontier_df,
    loss_col       = "frob_H",
    constraint_col = "var_H"
  )
  best_row   <- knee_res$best_row
  best_trial <- best_row$trial_index

  # ---- Plot ----------------------------------------------------------------
  plot_res <- .plot_pareto_kraljic(
    data_all             = trials_df,
    data_frontier        = frontier_df,
    bestPareto           = best_row,
    constraint_threshold = constraint_threshold
  )

  # ---- Output --------------------------------------------------------------
  list(
    hyperparameters = list(
      lambda = best_row$lambda_,
      gamma  = best_row$gamma
    ),
    metrics = list(
      loss       = best_row$loss,
      constraint = best_row$constraint
    ),
    trials = trials_df,
    W      = W[[best_trial]],
    H      = H[[best_trial]],
    plot   = plot_res
  )
}
