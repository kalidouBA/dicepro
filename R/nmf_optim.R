# =============================================================================
# NMF_optim.R
#
# Public API  : nmf_lbfgsb()
# Private fns : .norm_frob(), .truncate2()
#
# Rcpp entry points used (defined in src/):
#   compute_obj_term1_eigen_fast()
#   compute_grad_eigen_fast()
# =============================================================================


# -----------------------------------------------------------------------------
# nmf_lbfgsb  [public]
# -----------------------------------------------------------------------------

#' NMF with L-BFGS-B optimization
#'
#' Performs non-negative matrix factorization (NMF) using L-BFGS-B
#' optimization. The objective combines a Gaussian log-likelihood with an
#' augmented Lagrangian penalty that enforces the sum-to-one constraint on
#' cell-type proportions.
#'
#' @param r_dataset  Named list with three elements:
#' \describe{
#'   \item{B}{Numeric matrix (N_gene × N_sample). Bulk expression to
#'     factorize.}
#'   \item{P_cb}{Numeric matrix (N_sample × N_cells_Type). Known cell-type
#'     proportion estimates from a supervised deconvolution step.}
#'   \item{W_cb}{Numeric matrix (N_gene × N_cells_Type). Reference signature
#'     matrix for known cell types.}
#' }
#' @param W_prime    Numeric scalar or matrix. Initial value(s) for the
#'   unknown cell-type column(s) in \eqn{W}. Default \code{0}.
#' @param p_prime    Numeric scalar or matrix. Initial value(s) for the
#'   unknown cell-type column(s) in \eqn{H}. Default \code{0}.
#' @param lambda_    Numeric scalar. Augmented Lagrangian multiplier
#'   initialization (default \code{10}).
#' @param gamma_par  Numeric scalar. Penalty coefficient for the sum-to-one
#'   constraint (default \code{100}). Reduce if the optimizer returns
#'   infinite values.
#' @param path2save  Character scalar. Directory where intermediate results
#'   may be saved (currently unused; reserved for future use).
#' @param N_unknownCT Positive integer. Number of latent cell types to
#'   estimate (default \code{1}).
#' @param con        Named list of control parameters forwarded to
#'   \code{lbfgsb3c::lbfgsb3c}. Default: \code{list(maxit = 3000)}.
#'
#' @return A named list on success:
#' \describe{
#'   \item{W}{Optimized signature matrix (N_gene × N_cells_Type).}
#'   \item{frob_W}{Frobenius norm of \code{W}.}
#'   \item{var_W}{Variance of the unknown column of \code{W}.}
#'   \item{H}{Optimized proportion data.frame (N_sample × N_cells_Type + 1),
#'     with a leading \code{Mixture} column.}
#'   \item{frob_H}{Frobenius norm of \code{H}.}
#'   \item{var_H}{Variance of the unknown column of \code{H}.}
#'   \item{loss}{Gaussian log-likelihood (objective terms 1 + 2).}
#'   \item{frobNorm}{Objective term 1 (data-fit component).}
#'   \item{constNorm}{Objective term 2 (log-normalization component).}
#'   \item{c1}{Penalty term 1 (linear Lagrangian).}
#'   \item{c2}{Penalty term 2 (quadratic penalty).}
#'   \item{objectiveValue}{Total objective value (\code{loss}).}
#'   \item{penalty}{Total penalty (\code{c1 + c2}).}
#'   \item{cvrge}{Convergence code from \code{lbfgsb3c}.}
#'   \item{constraint}{Absolute constraint violation: \eqn{|1 - \sum h_i|}.}
#' }
#' Returns \code{NULL} when the optimizer fails or the solution is degenerate.
#'
#' @export
nmf_lbfgsb <- function(r_dataset,
                       W_prime     = 0,
                       p_prime     = 0,
                       lambda_     = 10,
                       gamma_par   = 100,
                       path2save   = "",
                       N_unknownCT = 1L,
                       con         = list(maxit = 3000L)) {

  B    <- as.matrix(r_dataset$B)
  p_cb <- as.matrix(r_dataset$P_cb)
  W_cb <- as.matrix(r_dataset$W_cb)

  N_gene      <- nrow(B)
  N_sample    <- ncol(B)
  N_cellsType <- ncol(W_cb) + N_unknownCT

  # ---- Initialise unknown columns ------------------------------------------
  p_cb_C <- matrix(p_prime, nrow = N_sample, ncol = N_unknownCT)
  W_cb_C <- matrix(W_prime, nrow = N_gene,   ncol = N_unknownCT)

  unknown_names          <- paste0("Unknown_", seq_len(N_unknownCT))
  colnames(p_cb_C)       <- unknown_names
  colnames(W_cb_C)       <- unknown_names

  # ---- Assemble H and W ----------------------------------------------------
  H <- as.matrix(cbind(p_cb, p_cb_C))
  H <- H / rowSums(H)

  W <- as.matrix(cbind(W_cb, W_cb_C))

  # ---- Initialise sigma via asinh-transformed residuals --------------------
  residuals_vec             <- as.vector(tcrossprod(W, H) - B)
  transformed_residuals_vec <- asinh(residuals_vec)
  sigma_par                 <- stats::sd(transformed_residuals_vec)

  theta      <- c(W_cb_C, H, sigma_par)
  lambda_par <- rep(lambda_, N_sample)

  # ---- Objective function --------------------------------------------------
  obj_fun <- function(theta) {
    W[, N_cellsType] <- theta[seq_len(N_gene)]
    H_loc <- matrix(
      theta[seq.int(N_gene + 1L, length(theta) - 1L)],
      nrow = N_sample, ncol = N_cellsType
    )
    sigma_loc <- theta[length(theta)]

    obj1 <- compute_obj_term1_eigen_fast(W, H_loc, B, sigma_loc)
    obj2 <- (N_sample * N_gene / 2) * log(2 * pi * sigma_loc^2)
    h_H  <- rowSums(H_loc) - 1
    obj3 <- lambda_par %*% h_H
    obj4 <- (gamma_par / 2) * sum(h_H^2)

    obj1 + obj2 + obj3 + obj4
  }

  # ---- Gradient function ---------------------------------------------------
  grad_obj_fun <- function(theta) {
    W[, N_cellsType] <- theta[seq_len(N_gene)]
    H_loc <- matrix(
      theta[seq.int(N_gene + 1L, length(theta) - 1L)],
      nrow = N_sample, ncol = N_cellsType
    )
    sigma_loc <- theta[length(theta)]

    compute_grad_eigen_fast(W, H_loc, B, sigma_loc, lambda_par, gamma_par)
  }

  # ---- Optimise ------------------------------------------------------------
  result <- try(
    lbfgsb3c::lbfgsb3c(
      par     = theta,
      fn      = obj_fun,
      gr      = grad_obj_fun,
      lower   = c(rep(0,    length(theta) - 1L), 1e-6),
      upper   = c(rep(Inf,  nrow(W_cb_C)),
                  rep(1,    N_sample * N_cellsType),
                  Inf),
      control = con
    ),
    silent = TRUE
  )

  # ---- Handle optimiser failure --------------------------------------------
  if (inherits(result, "try-error")) {
    message("Error from lbfgsb3c:\n", result[1L], "\n")
    warning(
      "Infinite values during optimization. Try lower values for 'gamma_par'.",
      call. = FALSE
    )
    return(NULL)
  }

  # ---- Extract solution ----------------------------------------------------
  theta <- result$par

  W_opt                 <- W
  W_opt[, N_cellsType]  <- theta[seq_len(N_gene)]

  H_opt <- matrix(
    theta[seq.int(N_gene + 1L, length(theta) - 1L)],
    nrow = N_sample, ncol = N_cellsType
  )
  dimnames(H_opt) <- dimnames(H)

  sigma_opt <- theta[length(theta)]

  # ---- Diagnostics ---------------------------------------------------------
  frob_H <- .norm_frob(H_opt)
  var_H  <- stats::var(as.vector(H_opt[, N_cellsType]))
  frob_W <- .norm_frob(W_opt)
  var_W  <- stats::var(as.vector(W_opt[, N_cellsType]))

  if (any(is.na(c(frob_H, var_H, frob_W, var_W))) ||
      any(c(var_H, var_W) == 0)) {
    return(NULL)
  }

  # ---- Constraint satisfaction ---------------------------------------------
  h_H         <- abs(rowSums(H_opt) - 1)
  constraints <- sum(h_H)

  # ---- Rebuild H as data.frame with Mixture column -------------------------
  H_opt[H_opt < 0] <- 0
  row_sums <- rowSums(H_opt)
  row_sums[row_sums == 0] <- 1
  H_opt <- H_opt / row_sums

  H_df <- as.data.frame(H_opt)
  H_df <- cbind(Mixture = rownames(H_df), H_df)


  # ---- Final objective decomposition ---------------------------------------
  obj1 <- compute_obj_term1_eigen_fast(W_opt, H_opt, B, sigma_opt)
  obj2 <- (N_sample * N_gene / 2) * log(2 * pi * sigma_opt^2)
  obj3 <- sum(lambda_par * h_H)
  obj4 <- (gamma_par / 2) * sum(h_H^2)

  list(
    W              = W_opt,
    frob_W         = frob_W,
    var_W          = var_W,
    H              = H_df,
    frob_H         = frob_H,
    var_H          = var_H,
    loss           = obj1 + obj2,
    frobNorm       = obj1,
    constNorm      = obj2,
    c1             = obj3,
    c2             = obj4,
    objectiveValue = obj1 + obj2,
    penalty        = obj3 + obj4,
    cvrge          = result$convergence,
    constraint     = abs(1 - constraints)
  )
}


# -----------------------------------------------------------------------------
# .truncate2  [private]
# -----------------------------------------------------------------------------

#' Truncate a numeric vector to two decimal places (no rounding)
#'
#' @param x Numeric vector.
#' @return Numeric vector truncated to two decimal places.
#' @keywords internal
#' @noRd
.truncate2 <- function(x) floor(x * 100) / 100


# -----------------------------------------------------------------------------
# .norm_frob  [private]
# -----------------------------------------------------------------------------

#' Frobenius norm of a matrix
#'
#' @param H Numeric matrix.
#' @return Non-negative numeric scalar.
#' @keywords internal
#' @noRd
.norm_frob <- function(H) sqrt(sum(as.matrix(H)^2))
