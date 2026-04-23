# =============================================================================
# Performance Metrics Functions
# =============================================================================


# -----------------------------------------------------------------------------
# row_norm_pos
# -----------------------------------------------------------------------------

#' Row-normalize a matrix with non-negative clamping
#'
#' Clamps negative values to 0, then divides each row by its sum so that
#' rows sum to 1. Rows whose sum equals 0 become all-\code{NaN} (degenerate
#' samples are flagged rather than silently zeroed).
#'
#' @param mat Numeric matrix.
#'
#' @return Numeric matrix of the same dimensions as \code{mat}, with
#'   each row summing to 1 (or all-\code{NaN} when the row sum is 0).
#'
#' @examples
#' m <- matrix(c(-1, 2, 3, 0, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
#' row_norm_pos(m)
#'
#' @export
row_norm_pos <- function(mat) {
  mat[mat < 0] <- 0
  rs <- rowSums(mat, na.rm = TRUE)
  rs[rs == 0] <- NA   # degenerate rows -> NaN, not silent zeros
  mat / rs
}


# -----------------------------------------------------------------------------
# samplewise_metrics
# -----------------------------------------------------------------------------

#' Sample-wise Pearson correlation and RMSE
#'
#' For each row (sample) of two matrices, computes the Pearson correlation
#' and RMSE between predicted and observed values, then returns their
#' cross-sample means.
#'
#' Rows where either \code{obs_mat} or \code{pred_mat} are constant
#' (zero standard deviation) yield \code{NA} for correlation but still
#' contribute to the RMSE mean.
#'
#' @param obs_mat  Numeric matrix of observations (rows = samples,
#'   columns = variables).
#' @param pred_mat Numeric matrix of predictions; must have the same
#'   dimensions as \code{obs_mat}.
#'
#' @return Named list with two elements:
#' \describe{
#'   \item{Correlation_mean}{Mean Pearson correlation across samples
#'     (ignoring \code{NA} rows).}
#'   \item{RMSE_mean}{Mean RMSE across samples.}
#' }
#'
#' @seealso \code{\link{full_metrics}}
#'
#' @export
samplewise_metrics <- function(obs_mat, pred_mat) {

  if (!identical(dim(obs_mat), dim(pred_mat)))
    stop("obs_mat and pred_mat must have identical dimensions.")

  N     <- nrow(obs_mat)
  cors  <- numeric(N)
  rmses <- numeric(N)

  for (i in seq_len(N)) {
    o <- obs_mat[i, ]
    p <- pred_mat[i, ]

    cors[i] <- if (sd(o, na.rm = TRUE) == 0 || sd(p, na.rm = TRUE) == 0) {
      NA_real_
    } else {
      cor(o, p, use = "pairwise.complete.obs")
    }

    rmses[i] <- sqrt(mean((o - p)^2, na.rm = TRUE))
  }

  list(
    Correlation_mean = mean(cors,  na.rm = TRUE),
    RMSE_mean        = mean(rmses, na.rm = TRUE)
  )
}


# -----------------------------------------------------------------------------
# full_metrics
# -----------------------------------------------------------------------------

#' Full agreement metrics via mixed-effects modeling
#'
#' Fits a two-way mixed-effects model (populations and subjects as random
#' effects) and derives ICC(3,1), Concordance correlation coefficient, and relative RMSE, complemented by
#' sample-wise correlation and RMSE from \code{\link{samplewise_metrics}}.
#'
#' @param x Numeric matrix of observed values (rows = subjects,
#'   columns = variables / populations).
#' @param y Numeric matrix of predicted values; same dimensions as \code{x}.
#'
#' @return Named list:
#' \describe{
#'   \item{Correlation_mean}{Mean sample-wise Pearson correlation.}
#'   \item{RMSE_mean}{Mean sample-wise RMSE.}
#'   \item{ICC3_adapted}{ICC(3,1) from the mixed model
#'     (\code{NA} for constant data).}
#'   \item{CCC_adapted}{Concordance correlation coefficient
#'     (\code{NA} for constant data).}
#'   \item{rRMSE_adapted}{Relative RMSE \eqn{\sqrt{\sigma^2_\varepsilon / V_T}}
#'     (\code{NA} for constant data).}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- matrix(rnorm(100), ncol = 5)
#' y <- x + matrix(rnorm(100, sd = 0.1), ncol = 5)
#' full_metrics(x, y)
#' }
#'
#' @export
full_metrics <- function(x, y) {

  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for full_metrics().")
  if (!identical(dim(x), dim(y)))
    stop("x and y must have identical dimensions.")

  cor_rmse <- samplewise_metrics(x, y)

  # Degenerate case: constant data -> skip mixed model
  if (sd(as.vector(x), na.rm = TRUE) == 0 ||
      sd(as.vector(y), na.rm = TRUE) == 0) {
    return(list(
      Correlation_mean = cor_rmse$Correlation_mean,
      RMSE_mean        = cor_rmse$RMSE_mean,
      ICC3_adapted     = NA_real_,
      CCC_adapted      = NA_real_,
      rRMSE_adapted    = NA_real_
    ))
  }

  n <- nrow(x)
  p <- ncol(x)

  df <- data.frame(
    value      = c(as.vector(x), as.vector(y)),
    subject    = factor(rep(seq_len(n), times = 2L * p)),
    population = factor(rep(rep(colnames(x), each = n), times = 2L)),
    method     = factor(rep(c("obs", "pred"), each = n * p))
  )

  model <- lme4::lmer(
    value ~ 1 + (1 | population / subject),
    data    = df,
    control = lme4::lmerControl(optimizer = "bobyqa")
  )

  vc             <- lme4::VarCorr(model)
  var_population <- vc[["population"]][1L, 1L]
  var_subject    <- vc[["subject:population"]][1L, 1L]
  var_residual   <- attr(vc, "sc")^2
  V_T            <- var_population + var_subject

  if (V_T == 0) {
    icc3  <- NA_real_
    ccc   <- NA_real_
    rrmse <- NA_real_
  } else {
    icc3  <- V_T / (V_T + var_residual)
    ccc   <- 2 * V_T / (2 * V_T + var_residual)
    rrmse <- sqrt(var_residual / V_T)
  }

  list(
    Correlation_mean = cor_rmse$Correlation_mean,
    RMSE_mean        = cor_rmse$RMSE_mean,
    ICC3_adapted     = icc3,
    CCC_adapted      = ccc,
    rRMSE_adapted    = rrmse
  )
}


# -----------------------------------------------------------------------------
# MakeTable1Tool
# -----------------------------------------------------------------------------

#' Build a performance metrics table for composition matrices
#'
#' Aligns predicted and observed matrices to their common populations,
#' row-normalizes both, then computes global and sample-wise agreement
#' metrics.
#'
#' @param pred_mat Numeric matrix of predicted compositions
#'   (rows = samples, columns = populations).
#' @param obs_mat  Numeric matrix of observed compositions; must share at
#'   least one column name with \code{pred_mat}.
#' @importFrom stats cor
#' @return A list with one element:
#' \describe{
#'   \item{Perf}{A one-row \code{data.frame} containing:
#'     \code{Correlation}, \code{RMSE}, \code{Correlation_mean},
#'     \code{RMSE_mean}, \code{NRMSE}, \code{ICC3}, \code{CCC}.
#'     All fields are \code{NA} when computation is not possible.}
#' }
#'
#' @details
#' Populations absent from \code{pred_mat} but present in \code{obs_mat}
#' (and vice versa) are silently dropped after the intersection step.
#' Both matrices are row-normalized via \code{\link{row_norm_pos}} before
#' any metric is computed.
#'
#' @seealso \code{\link{full_metrics}}, \code{\link{row_norm_pos}}
#'
#' @export
MakeTable1Tool <- function(pred_mat, obs_mat) {

  na_row <- data.frame(
    Correlation      = NA_real_,
    RMSE             = NA_real_,
    Correlation_mean = NA_real_,
    RMSE_mean        = NA_real_,
    NRMSE            = NA_real_,
    ICC3             = NA_real_,
    CCC              = NA_real_
  )

  pred_mat <- as.matrix(pred_mat)
  obs_mat  <- as.matrix(obs_mat)

  # ---- Align to common populations -------------------------------------------
  common <- intersect(colnames(pred_mat), colnames(obs_mat))
  if (length(common) == 0L)
    stop("No common populations between pred_mat and obs_mat.")

  pred_mat <- pred_mat[, common, drop = FALSE]
  obs_mat  <- obs_mat[,  common, drop = FALSE]

  # ---- Row-normalize ----------------------------------------------------------
  pred_mat <- row_norm_pos(pred_mat)
  obs_mat  <- row_norm_pos(obs_mat)

  obs_vec  <- as.vector(t(obs_mat))
  pred_vec <- as.vector(t(pred_mat))

  if (sd(pred_vec, na.rm = TRUE) == 0) {
    message("Skipping: constant prediction vector.")
    return(list(Perf = na_row))
  }

  # ---- Global metrics --------------------------------------------------------
  na_row$Correlation <- cor(obs_vec, pred_vec, use = "pairwise.complete.obs")
  na_row$RMSE        <- sqrt(mean((obs_vec - pred_vec)^2, na.rm = TRUE))

  # ---- Mixed-model metrics ---------------------------------------------------
  m <- full_metrics(x = obs_mat, y = pred_mat)

  na_row$Correlation_mean <- m$Correlation_mean
  na_row$RMSE_mean        <- m$RMSE_mean
  na_row$NRMSE            <- m$rRMSE_adapted
  na_row$ICC3             <- m$ICC3_adapted
  na_row$CCC              <- m$CCC_adapted

  list(Perf = na_row)
}



# -----------------------------------------------------------------------------
# nrmse
# -----------------------------------------------------------------------------

#' Normalized Root Mean Square Error (NRMSE)
#'
#' Computes RMSE between back-transformed \code{pred} and \code{obs} vectors,
#' then normalizes by a summary statistic of the observations.
#'
#' @param pred           Numeric vector of predicted values (on the
#'   transformed scale).
#' @param obs            Numeric vector of observed values (same scale as
#'   \code{pred}; same length).
#' @param method         Normalization denominator: \code{"sd"} (default),
#'   \code{"mean"}, \code{"maxmin"}, or \code{"iq"} (interquartile range).
#' @param transformation Back-transformation applied before computing RMSE.
#'   One of \code{"none"} (default), \code{"sqrt"}, \code{"4thrt"},
#'   \code{"log"}, \code{"log10"}, \code{"log2"}, \code{"log1p"},
#'   \code{"arcsine"}, or \code{"other"}.
#' @param trans_function R expression string used when
#'   \code{transformation = "other"} (the variable inside the expression
#'   must be named \code{x}).
#'
#' @return Non-negative numeric scalar, or \code{NA_real_} (with a message)
#'   when NRMSE is undefined (fewer than 2 complete observations, or zero
#'   denominator).
#'
#' @examples
#' set.seed(1)
#' obs  <- rnorm(100)
#' pred <- obs + rnorm(100, sd = 0.1)
#' nrmse(pred, obs)
#' nrmse(pred, obs, method = "maxmin")
#'
#' @importFrom stats sd quantile
#' @export
nrmse <- function(pred, obs,
                  method         = "sd",
                  transformation = "none",
                  trans_function = "none") {

  # ---- Input validation ------------------------------------------------------
  if (length(pred) != length(obs))
    stop(sprintf(
      "pred and obs must have the same length (obs: %d, pred: %d).",
      length(obs), length(pred)
    ))

  method <- match.arg(method, c("sd", "mean", "maxmin", "iq"))

  valid_trans <- c("none", "sqrt", "4thrt", "log", "log10",
                   "log2", "log1p", "arcsine", "other")
  transformation <- match.arg(transformation, valid_trans)

  if (transformation == "other" && identical(trans_function, "none"))
    stop("trans_function must be specified when transformation = 'other'.")

  # ---- Back-transformation ---------------------------------------------------
  back <- switch(
    transformation,
    sqrt    = function(v) v^2,
    `4thrt` = function(v) v^4,
    log     = exp,
    log10   = function(v) 10^v,
    log2    = function(v) 2^v,
    log1p   = expm1,
    arcsine = function(v) sin(v)^2,
    other   = function(v) { x <- v; eval(parse(text = trans_function)) },
    none    = identity
  )

  obs  <- back(obs)
  pred <- back(pred)

  # ---- Remove incomplete pairs -----------------------------------------------
  keep <- !is.na(obs) & !is.na(pred)
  obs  <- obs[keep]
  pred <- pred[keep]

  if (length(obs) < 2L) {
    message("Not enough observations to calculate NRMSE -- returning NA.")
    return(NA_real_)
  }

  # ---- RMSE ------------------------------------------------------------------
  rmse <- sqrt(mean((obs - pred)^2))

  # ---- Normalization denominator ---------------------------------------------
  denom <- switch(
    method,
    sd     = sd(obs),
    mean   = mean(obs),
    maxmin = diff(range(obs)),
    iq     = as.numeric(diff(quantile(obs, c(0.25, 0.75))))
  )

  if (is.na(denom) || denom == 0) {
    message(sprintf(
      "Denominator for method '%s' is zero or NA -- returning NA.", method
    ))
    return(NA_real_)
  }

  abs(rmse / denom)
}


