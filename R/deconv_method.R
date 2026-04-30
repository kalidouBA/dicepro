# =============================================================================
# deconv_method.R
#
# Public API  : running_method()
# Private fns : .run_cibersort_core(), .run_cibersort_perm(),
#               .cibersort_svr()
#
# Supported methods:
#   "CSx"        - CIBERSORTx via Docker (requires credentials)
#   "CS"  - CIBERSORT nu-SVR (built-in, no credentials needed)
#   "DCQ"        - ComICS::dcq
#   "FARDEEP"    - FARDEEP::fardeep
# =============================================================================


# =============================================================================
# CIBERSORT internal implementation
# Adapted from Newman et al. (2015), Stanford University.
# Original script: CIBERSORT R v1.04 (Aaron M. Newman, amnewman@stanford.edu)
# Licence: http://cibersort.stanford.edu/CIBERSORT_License.txt
# =============================================================================


# -----------------------------------------------------------------------------
# .cibersort_svr  [private]
# -----------------------------------------------------------------------------

#' Fit nu-SVR models and return the best-fit weights
#'
#' Tries three values of nu (0.25, 0.50, 0.75) in parallel, selects the
#' model with the lowest RMSE, and returns normalized cell-type weights.
#'
#' @param X          Numeric matrix. Standardized signature matrix
#'   (genes x cell types).
#' @param y          Numeric vector. Standardized mixture sample (length =
#'   number of genes).
#' @param absolute   Logical. Run in absolute mode.
#' @param abs_method Character. \code{"sig.score"} or \code{"no.sumto1"}.
#'
#' @return Named list: \code{w} (weights), \code{mix_rmse}, \code{mix_r}.
#' @keywords internal
#' @noRd
.cibersort_svr <- function(X, y, absolute, abs_method) {

  nu_vals <- c(0.25, 0.50, 0.75)

  # Fit one model per nu value - single-core on Windows, parallel elsewhere
  n_cores <- if (.Platform$OS.type == "windows") 1L else length(nu_vals)

  models <- parallel::mclapply(nu_vals, function(nu) {
    e1071::svm(X, y,
               type   = "nu-regression",
               kernel = "linear",
               nu     = nu,
               scale  = FALSE)
  }, mc.cores = n_cores)

  # Evaluate each model
  rmses <- vapply(models, function(m) {
    w_raw         <- t(m$coefs) %*% m$SV
    w_raw[w_raw < 0] <- 0
    w_norm        <- w_raw / sum(w_raw)
    fitted        <- rowSums(sweep(X, 2L, w_norm, `*`))
    sqrt(mean((fitted - y)^2))
  }, numeric(1L))

  corrs <- vapply(models, function(m) {
    w_raw         <- t(m$coefs) %*% m$SV
    w_raw[w_raw < 0] <- 0
    w_norm        <- w_raw / sum(w_raw)
    fitted        <- rowSums(sweep(X, 2L, w_norm, `*`))
    stats::cor(fitted, y)
  }, numeric(1L))

  best    <- which.min(rmses)
  model   <- models[[best]]

  # Extract and normalise final weights
  q         <- t(model$coefs) %*% model$SV
  q[q < 0]  <- 0

  w <- if (!absolute || abs_method == "sig.score") {
    q / sum(q)          # relative fractions
  } else {
    q                   # absolute scores (no.sumto1)
  }

  list(w = w, mix_rmse = rmses[best], mix_r = corrs[best])
}


# -----------------------------------------------------------------------------
# .run_cibersort_perm  [private]
# -----------------------------------------------------------------------------

#' Build empirical null distribution of correlations via permutation
#'
#' @param perm       Integer. Number of permutations.
#' @param X          Numeric matrix. Standardised signature matrix.
#' @param Y          Numeric matrix. Mixture matrix (genes x samples).
#' @param absolute   Logical.
#' @param abs_method Character.
#'
#' @return Numeric vector of permuted correlations (sorted ascending).
#' @keywords internal
#' @noRd
.run_cibersort_perm <- function(perm, X, Y, absolute, abs_method) {

  all_y <- as.vector(Y)

  null_r <- vapply(seq_len(perm), function(.) {
    yr <- sample(all_y, nrow(X), replace = FALSE)
    yr <- (yr - mean(yr)) / stats::sd(yr)
    .cibersort_svr(X, yr, absolute, abs_method)$mix_r
  }, numeric(1L))

  sort(null_r)
}


# -----------------------------------------------------------------------------
# .run_cibersort_core  [private]
# -----------------------------------------------------------------------------

#' Run CIBERSORT deconvolution on intersected matrices
#'
#' Implements the full CIBERSORT pipeline: optional quantile normalisation,
#' gene intersection, signature standardisation, per-sample nu-SVR, optional
#' permutation p-values, and optional absolute scoring.
#'
#' @param sig_matrix   Numeric matrix. Reference signature (genes x cell types).
#' @param mixture      Numeric matrix. Bulk expression (genes x samples).
#' @param perm         Non-negative integer. Permutations for p-values
#'   (default \code{0L} = no p-values).
#' @param QN           Logical. Quantile-normalise the mixture (default
#'   \code{TRUE}).
#' @param absolute     Logical. Return absolute scores instead of fractions
#'   (default \code{FALSE}).
#' @param abs_method   Character. \code{"sig.score"} (default) or
#'   \code{"no.sumto1"}.
#'
#' @return Numeric matrix (samples x cell types + metadata columns).
#'   Columns: one per cell type, then \code{P.value}, \code{Correlation},
#'   \code{RMSE} (and \code{Absolute.score} when \code{absolute = TRUE}).
#'
#' @keywords internal
#' @noRd
.run_cibersort_core <- function(sig_matrix,
                                mixture,
                                perm       = 0L,
                                QN         = TRUE,
                                absolute   = FALSE,
                                abs_method = "sig.score") {

  if (absolute && !abs_method %in% c("sig.score", "no.sumto1"))
    stop("abs_method must be 'sig.score' or 'no.sumto1'.", call. = FALSE)

  X <- data.matrix(sig_matrix)
  Y <- data.matrix(mixture)

  # ---- Deduplicate rownames in Y ------------------------------------------
  dups <- nrow(Y) - length(unique(rownames(Y)))
  if (dups > 0L) {
    warning(sprintf("%d duplicated gene symbol(s) in mixture - made unique.",
                    dups), call. = FALSE)
    rownames(Y) <- make.names(rownames(Y), unique = TRUE)
  }

  # ---- Sort rows -----------------------------------------------------------
  X <- X[order(rownames(X)), , drop = FALSE]
  Y <- Y[order(rownames(Y)), , drop = FALSE]

  # ---- Anti-log if mixture looks log-transformed ---------------------------
  if (max(Y, na.rm = TRUE) < 50) Y <- 2^Y

  # ---- Quantile normalisation ----------------------------------------------
  if (QN) {
    cn <- colnames(Y); rn <- rownames(Y)
    Y  <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- cn; rownames(Y) <- rn
  }

  Yorig   <- Y
  Ymedian <- max(stats::median(Yorig), 1)

  # ---- Gene intersection ---------------------------------------------------
  common <- intersect(rownames(X), rownames(Y))
  if (length(common) == 0L)
    stop("No common genes between signature matrix and mixture.", call. = FALSE)
  X <- X[common, , drop = FALSE]
  Y <- Y[common, , drop = FALSE]

  # ---- Standardise signature matrix ----------------------------------------
  X <- (X - mean(X)) / stats::sd(as.vector(X))

  # ---- Empirical null distribution (permutations) -------------------------
  null_dist <- if (perm > 0L) {
    .run_cibersort_perm(perm, X, Y, absolute, abs_method)
  } else {
    NULL
  }

  # ---- Per-sample SVR ------------------------------------------------------
  n_samples <- ncol(Y)

  results <- lapply(seq_len(n_samples), function(i) {

    y <- Y[, i]
    y <- (y - mean(y)) / stats::sd(y)

    fit      <- .cibersort_svr(X, y, absolute, abs_method)
    w        <- fit$w
    mix_r    <- fit$mix_r
    mix_rmse <- fit$mix_rmse

    # Absolute sig.score scaling
    if (absolute && abs_method == "sig.score") {
      w <- w * stats::median(Yorig[, i]) / Ymedian
    }

    # P-value from null distribution
    pval <- if (!is.null(null_dist)) {
      1 - which.min(abs(null_dist - mix_r)) / length(null_dist)
    } else {
      NA_real_
    }

    row <- c(as.vector(w), pval, mix_r, mix_rmse)
    if (absolute) row <- c(row, sum(w))
    row
  })

  # ---- Assemble output matrix ----------------------------------------------
  out <- do.call(rbind, results)
  rownames(out) <- colnames(Y)

  ct_names <- colnames(X)
  col_names <- c(ct_names, "P.value", "Correlation", "RMSE")
  if (absolute)
    col_names <- c(col_names,
                   paste0("Absolute.score.", abs_method))
  colnames(out) <- col_names

  out
}


# =============================================================================
# running_method  [public]
# =============================================================================

#' Run cell-type deconvolution
#'
#' Performs cell-type deconvolution on bulk RNA-seq data using one of four
#' supported methods. All methods receive gene-intersected, numeric matrices
#' and return a proportion matrix normalised to sum to 1 per sample
#' (columns).
#'
#' @param bulk             Numeric matrix (genes x samples). Bulk RNA-seq
#'   expression data. Row names must be gene symbols.
#' @param reference        Numeric matrix (genes x cell types). Reference
#'   cell-type signature profiles. Row names must be gene symbols.
#' @param methodDeconv     Character scalar. Deconvolution method:
#' \describe{
#'   \item{\code{"CSx"}}{CIBERSORTx - requires Docker and credentials.}
#'   \item{\code{"CS"}}{Built-in CIBERSORT nu-SVR implementation.
#'     Requires \pkg{e1071}, \pkg{parallel}, and \pkg{preprocessCore}.
#'     No external account needed.}
#'   \item{\code{"DCQ"}}{Requires \pkg{ComICS}.}
#'   \item{\code{"FARDEEP"}}{Requires \pkg{FARDEEP}.}
#' }
#' @param cibersortx_email Character. CIBERSORTx account email (required
#'   when \code{methodDeconv = "CSx"}).
#' @param cibersortx_token Character. CIBERSORTx account token (required
#'   when \code{methodDeconv = "CSx"}).
#' @param cibersort_perm   Non-negative integer. Number of permutations for
#'   CIBERSORT p-values (default \code{0L} = no p-values). Only used when
#'   \code{methodDeconv = "CS"}.
#' @param cibersort_QN     Logical. Apply quantile normalisation in
#'   CIBERSORT (default \code{TRUE}). Only used when
#'   \code{methodDeconv = "CS"}.
#'
#' @return Numeric matrix (cell types x samples) of estimated cell-type
#'   proportions. Each column sums to 1.
#'
#' @details
#' Gene intersection is performed automatically: only genes present in both
#' \code{bulk} and \code{reference} are used. Both matrices are coerced to
#' numeric before processing.
#'
#' For \code{"CS"}, the built-in implementation closely follows the
#' original Newman et al. (2015) algorithm: nu-SVR with three candidate nu
#' values (0.25, 0.50, 0.75), optional quantile normalisation, and optional
#' permutation-based p-values. Packages \pkg{e1071}, \pkg{parallel}, and
#' \pkg{preprocessCore} must be installed.
#'
#' @seealso \code{\link{run_CSx}}
#'
#' @references
#' Newman AM et al. (2015). Robust enumeration of cell subsets from tissue
#' expression profiles. \emph{Nature Methods}, 12(5), 453-457.
#' \doi{10.1038/nmeth.3337}
#'
#' @export
running_method <- function(bulk,
                           reference,
                           methodDeconv     = "CSx",
                           cibersortx_email = NULL,
                           cibersortx_token = NULL,
                           cibersort_perm   = 0L,
                           cibersort_QN     = TRUE) {

  valid_methods <- c("CSx", "CS", "DCQ", "FARDEEP")

  methodDeconv <- tryCatch(
    match.arg(methodDeconv, valid_methods),
    error = function(e) stop(sprintf(
      "Invalid 'methodDeconv': '%s'. Must be one of: %s.",
      methodDeconv, paste(valid_methods, collapse = ", ")
    ), call. = FALSE)
  )

  # ---- Coerce to numeric matrices preserving dimnames ---------------------
  .to_numeric_matrix <- function(mat) {
    dn  <- dimnames(mat)
    mat <- apply(mat, 2L, as.numeric)
    dimnames(mat) <- dn
    mat
  }
  reference <- .to_numeric_matrix(reference)
  bulk      <- .to_numeric_matrix(bulk)

  # ---- Gene intersection --------------------------------------------------
  common <- intersect(rownames(bulk), rownames(reference))
  if (length(common) == 0L)
    stop("No common genes between 'bulk' and 'reference'.", call. = FALSE)

  bulk_int <- bulk[common, , drop = FALSE]
  ref_int  <- reference[common, , drop = FALSE]

  # ---- Dispatch -----------------------------------------------------------
  out_Dec <- switch(
    methodDeconv,

    CSx = {
      if (is.null(cibersortx_email) || is.null(cibersortx_token))
        stop("CIBERSORTx credentials required for methodDeconv = 'CSx'.",
             call. = FALSE)
      run_CSx(bulk, reference, cibersortx_email, cibersortx_token)
    },

    CS = {
      # Check required packages
      for (pkg in c("e1071", "parallel", "preprocessCore")) {
        if (!requireNamespace(pkg, quietly = TRUE))
          stop(sprintf(
            "Package '%s' is required for methodDeconv = 'CS'.\n",
            pkg
          ), call. = FALSE)
      }
      raw <- .run_cibersort_core(
        sig_matrix = ref_int,
        mixture    = bulk_int,
        perm       = cibersort_perm,
        QN         = cibersort_QN,
        absolute   = FALSE,
        abs_method = "sig.score"
      )
      # Keep only cell-type columns (drop P.value, Correlation, RMSE)
      ct_cols <- colnames(ref_int)
      t(raw[, ct_cols, drop = FALSE])
    },

    DCQ = {
      if (!requireNamespace("ComICS", quietly = TRUE))
        stop("Package 'ComICS' is required for methodDeconv = 'DCQ'.",
             call. = FALSE)
      raw <- ComICS::dcq(
        reference_data    = as.data.frame(ref_int),
        mix_data          = as.data.frame(bulk_int),
        marker_set        = as.data.frame(rownames(ref_int)),
        alpha_used        = 0.05,
        lambda_min        = 0.2,
        number_of_repeats = 30L
      )$average
      raw <- t(raw)
      raw[raw < 0] <- 0
      raw
    },

    FARDEEP = {
      if (!requireNamespace("FARDEEP", quietly = TRUE))
        stop("Package 'FARDEEP' is required for methodDeconv = 'FARDEEP'.",
             call. = FALSE)
      suppressWarnings(
        t(FARDEEP::fardeep(
          X         = reference,
          Y         = bulk,
          nn        = TRUE,
          intercept = TRUE,
          permn     = 100L,
          QN        = FALSE
        )$abs.beta)
      )
    }
  )

  # ---- Normalise columns to sum to 1 --------------------------------------
  col_sums <- colSums(out_Dec)
  col_sums[col_sums == 0] <- 1
  out_Dec <- sweep(out_Dec, 2L, col_sums, FUN = "/")

  out_Dec
}
