# =============================================================================
# hypersearch.R
#
# Public API  : contains_nan_or_inf(), nmf_lbfgsb_hyperOpt(),
#               objective_opt(), research_hyperOpt()
# Private fns : .clean_nmf_matrix(), objective_wrapper(),
#               .sample_from_space(), .custom_progress_bar(),
#               .tpe_sample(), .kde_log_density(),
#               .parse_config(), .parse_hyperopt_searchspace(),
#               .get_report_path()
# =============================================================================


# -----------------------------------------------------------------------------
# contains_nan_or_inf  [public]
# -----------------------------------------------------------------------------

#' Check if a value contains NaN or Inf
#'
#' @param value Numeric, vector, matrix, or data frame to check.
#' @return Logical. TRUE if NaN or Inf is present, FALSE otherwise.
#' @export
contains_nan_or_inf <- function(value) {
  if (is.numeric(value)) {
    return(any(is.nan(value) | is.infinite(value)))
  } else if (is.matrix(value) || is.data.frame(value)) {
    return(any(is.nan(as.matrix(value)) | is.infinite(as.matrix(value))))
  }
  return(FALSE)
}


# -----------------------------------------------------------------------------
# .clean_nmf_matrix  [private]
# -----------------------------------------------------------------------------

#' Clean an H or W matrix returned by nmf_lbfgsb
#'
#' Strips any character columns (e.g. row-name columns stored as data),
#' promotes them to \code{rownames}, and coerces all remaining entries to
#' numeric.
#'
#' @param mat Matrix or data frame returned by \code{nmf_lbfgsb}.
#' @return A numeric matrix with correct \code{rownames} and \code{colnames}.
#' @keywords internal
#' @noRd
.clean_nmf_matrix <- function(mat) {
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (!is.matrix(mat)) return(mat)

  char_cols <- which(
    apply(mat, 2L, function(x)
      !suppressWarnings(all(!is.na(as.numeric(x))))
    )
  )

  if (length(char_cols) > 0L) {
    rn  <- mat[, char_cols[1L]]
    mat <- mat[, -char_cols, drop = FALSE]
    rownames(mat) <- rn
  }

  rn  <- rownames(mat)
  cn  <- colnames(mat)
  mat <- matrix(
    as.numeric(mat),
    nrow     = nrow(mat),
    ncol     = ncol(mat),
    dimnames = list(rn, cn)
  )
  mat
}


# -----------------------------------------------------------------------------
# nmf_lbfgsb_hyperOpt  [public]
# -----------------------------------------------------------------------------

#' NMF L-BFGS-B wrapper for hyperparameter optimisation
#'
#' Thin wrapper around \code{\link{nmf_lbfgsb}} that normalises the
#' \code{p_prime} argument and cleans the returned matrices.
#'
#' @param dataset   List containing matrices \code{B}, \code{W}, and \code{P}.
#' @param W_prime   Optional numeric matrix. Initial \eqn{W'}.
#' @param p_prime   Optional numeric matrix or scalar. Initial \eqn{P'}.
#'   If \code{NULL}, defaults to a matrix of \code{0.1}.
#' @param lambda_   Numeric. Regularisation parameter \eqn{\lambda}.
#' @param gamma_par Numeric. Regularisation parameter \eqn{\gamma}.
#' @param path2save Character. Path passed to \code{nmf_lbfgsb}.
#'
#' @return List output from \code{\link{nmf_lbfgsb}} with \code{H} and
#'   \code{W} cleaned by \code{.clean_nmf_matrix}.
#'
#' @export
nmf_lbfgsb_hyperOpt <- function(dataset,
                                W_prime   = NULL,
                                p_prime   = NULL,
                                lambda_   = 10,
                                gamma_par = 100,
                                path2save = "") {

  B    <- as.data.frame(dataset$B)
  W_cb <- as.data.frame(dataset$W)
  P_cb <- as.data.frame(dataset$P)

  # --- Normalise p_prime to a matrix ----------------------------------------
  if (is.null(p_prime)) {
    N_sample    <- ncol(B)
    N_unknownCT <- 1L
    p_prime     <- matrix(0.1, nrow = N_sample, ncol = N_unknownCT)

  } else if (is.numeric(p_prime) && is.vector(p_prime)) {
    N_sample    <- ncol(B)
    N_unknownCT <- 1L
    p_prime     <- matrix(p_prime, nrow = N_sample, ncol = N_unknownCT)

  } else if (is.matrix(p_prime)) {
    N_sample    <- nrow(p_prime)
    N_unknownCT <- ncol(p_prime)

  } else {
    stop("'p_prime' must be NULL, a numeric vector, or a matrix.", call. = FALSE)
  }

  # --- Ensure P_cb is (samples x cell_types) --------------------------------
  if (nrow(P_cb) != ncol(B)) P_cb <- t(P_cb)

  r_dataset <- list(B = B, W_cb = W_cb, P_cb = P_cb)

  result <- nmf_lbfgsb(
    r_dataset   = r_dataset,
    W_prime     = W_prime,
    p_prime     = p_prime,
    lambda_     = lambda_,
    gamma_par   = gamma_par,
    path2save   = path2save,
    N_unknownCT = N_unknownCT
  )

  if (!is.null(result$H)) result$H <- .clean_nmf_matrix(result$H)
  if (!is.null(result$W)) result$W <- .clean_nmf_matrix(result$W)

  return(result)
}


# -----------------------------------------------------------------------------
# objective_opt  [public]
# -----------------------------------------------------------------------------

#' Objective function for hyperparameter optimisation
#'
#' Runs one NMF trial for a given \eqn{(\lambda, \gamma, p')} configuration
#' and returns a structured result list, or \code{NULL} when the trial
#' produces invalid values.
#'
#' @param dataset      List with matrices \code{B}, \code{W}, and \code{P}.
#' @param config       List. Configuration object (needs \code{$exp} for the
#'   output directory).
#' @param lambda_      Numeric. Regularisation parameter \eqn{\lambda}.
#' @param gamma_factor Numeric or \code{NULL}. When not \code{NULL},
#'   \eqn{\gamma} is derived as \code{lambda_ * gamma_factor}.
#' @param gamma        Numeric. \eqn{\gamma} used directly when
#'   \code{gamma_factor} is \code{NULL}.
#' @param p_prime      Numeric matrix (\code{nrow = n_samples},
#'   \code{ncol = n_unknown_ct}).
#' @param W_prime      Numeric matrix or scalar. Initial \eqn{W'}.
#'
#' @return A named list with elements \code{loss}, \code{constraint},
#'   \code{status}, \code{current_params}, \code{W}, \code{H}, \code{cvrge};
#'   or \code{NULL} when the trial is invalid.
#'
#' @export
objective_opt <- function(dataset,
                          config       = list(),
                          lambda_      = NULL,
                          gamma_factor = NULL,
                          gamma        = NULL,
                          p_prime      = NULL,
                          W_prime      = 0) {

  exp_dir <- if (!is.null(config$exp)) config$exp else "."

  if (!is.null(gamma_factor)) {
    gamma <- lambda_ * gamma_factor
  }

  result <- nmf_lbfgsb_hyperOpt(
    dataset   = dataset,
    W_prime   = W_prime,
    p_prime   = p_prime,
    lambda_   = lambda_,
    gamma_par = gamma,
    path2save = exp_dir
  )

  # --- Coerce result components to plain scalars/matrices -------------------
  result_dict <- lapply(names(result), function(nm) {
    x <- result[[nm]]
    if (nm %in% c("H", "W", "p_prime_estm")) {
      if (is.data.frame(x)) as.matrix(x) else x
    } else if (is.numeric(x) && (is.null(dim(x)) || length(dim(x)) == 0L)) {
      x
    } else if (is.matrix(x) || is.data.frame(x)) {
      num_cols <- vapply(as.data.frame(x), is.numeric, logical(1L))
      as.numeric(as.matrix(as.data.frame(x)[, num_cols, drop = FALSE]))
    } else {
      x
    }
  })
  names(result_dict) <- names(result)

  # Helper: safely extract first element as numeric scalar
  .sc <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    as.numeric(x[[1L]])
  }

  # --- Guard: discard trials with non-finite values -------------------------
  if (contains_nan_or_inf(.sc(result_dict$objectiveValue)) ||
      contains_nan_or_inf(.sc(result_dict$constraint))     ||
      any(is.na(c(
        .sc(result_dict$frob_H), .sc(result_dict$var_H),
        .sc(result_dict$frob_W), .sc(result_dict$var_W)
      )))) {
    return(NULL)
  }

  list(
    loss       = .sc(result_dict$objectiveValue),
    constraint = .sc(result_dict$constraint),
    status     = "OK",
    current_params = list(
      gamma          = gamma,
      lambda_        = lambda_,
      p_prime        = p_prime[1L, 1L],
      frob_H         = .sc(result_dict$frob_H),
      var_H          = .sc(result_dict$var_H),
      frob_W         = .sc(result_dict$frob_W),
      var_W          = .sc(result_dict$var_W),
      frobNorm       = .sc(result_dict$frobNorm),
      constNorm      = .sc(result_dict$constNorm),
      c1             = .sc(result_dict$c1),
      c2             = .sc(result_dict$c2),
      objectiveValue = .sc(result_dict$objectiveValue),
      penalty        = .sc(result_dict$penalty)
    ),
    W     = result_dict$W,
    H     = result_dict$H,
    cvrge = result_dict$cvrge
  )
}


# -----------------------------------------------------------------------------
# objective_wrapper  [private]
# -----------------------------------------------------------------------------

#' Objective wrapper for one hyperparameter trial
#'
#' Calls \code{objective_opt()} inside a \code{tryCatch}, suppresses console
#' output, records timing, and returns a tidy \code{list(df, W, H)} or
#' \code{NULL} on failure.
#'
#' @param objective_opt Function. The objective to call (passed explicitly to
#'   avoid relying on lexical scoping across files).
#' @param dataset       List with \code{$B}, \code{$W}, \code{$P}.
#' @param config        List. Parsed configuration object.
#' @param params        List. Must contain \code{lambda_}, \code{gamma},
#'   \code{p_prime}; optionally \code{gamma_factor}.
#' @param W_prime       Matrix or scalar. Optional initial \eqn{W} for NMF.
#'
#' @return \code{list(df, W, H)} on success; \code{NULL} on failure.
#'
#' @importFrom utils capture.output
#' @keywords internal
#' @noRd
objective_wrapper <- function(objective_opt, dataset, config, params,
                              W_prime = NULL) {
  tryCatch(
    {
      start_time <- Sys.time()

      n_samples    <- ncol(dataset$B)
      n_unknown_ct <- if (
        !is.null(W_prime) && is.matrix(W_prime) && ncol(W_prime) > 0L
      ) ncol(W_prime) else 1L

      if (!is.matrix(params$p_prime)) {
        params$p_prime <- matrix(
          params$p_prime,
          nrow = n_samples,
          ncol = n_unknown_ct
        )
      }

      utils::capture.output(
        invisible(
          returned_dict <- objective_opt(
            dataset      = dataset,
            config       = config,
            lambda_      = params$lambda_,
            gamma_factor = params$gamma_factor,
            gamma        = params$gamma,
            p_prime      = params$p_prime,
            W_prime      = W_prime
          )
        ),
        file = nullfile()
      )

      if (is.null(returned_dict) ||
          (!is.null(returned_dict$status) && returned_dict$status != "OK")) {
        return(NULL)
      }

      duration <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      row <- list(
        lambda_      = params$lambda_,
        gamma        = params$gamma,
        gamma_factor = params$gamma_factor %||% NA_real_,
        p_prime      = params$p_prime[1L, 1L],
        loss         = returned_dict$loss,
        constraint   = returned_dict$constraint,
        status       = returned_dict$status,
        duration     = duration
      )

      if (!is.null(returned_dict$current_params)) {
        cp <- returned_dict$current_params
        cp <- cp[!names(cp) %in% c("gamma", "lambda_", "p_prime")]
        cp <- lapply(cp, function(v) {
          if (is.null(v) || length(v) == 0L) return(NA_real_)
          if (length(v) == 1L) return(v)
          v[[1L]]
        })
        row <- c(row, cp)
      }

      list(
        df = as.data.frame(row, check.names = FALSE),
        W  = returned_dict$W,
        H  = returned_dict$H
      )
    },
    error = function(e) {
      message("objective_wrapper error: ", conditionMessage(e))
      NULL
    }
  )
}


# -----------------------------------------------------------------------------
# .sample_from_space  [private]
# -----------------------------------------------------------------------------

#' Sample one configuration from a hyperparameter search space
#'
#' Each element of \code{space} is a named list with a \code{type} field and
#' the corresponding bounds/parameters. Atomic vectors of the form
#' \code{c(type, low, high)} are also accepted and parsed automatically.
#'
#' When both \code{lambda_factor} and \code{gamma} are present in the sampled
#' configuration, \code{lambda_} is derived as
#' \code{lambda_ = gamma * lambda_factor}.
#'
#' When \code{config$gamma_ratio_min} is set (mode \code{"all_gamma_dominant"}),
#' candidates where \code{gamma <= gamma_ratio_min * lambda_} are rejected and
#' resampled. After \code{max_tries} attempts the last candidate is returned
#' with a warning.
#'
#' @param space         Named list. Each element describes one hyperparameter.
#' @param gamma_ratio_min Numeric or \code{NULL}. Minimum ratio
#'   \code{gamma / lambda_} required. \code{NULL} disables the constraint
#'   (equivalent to the plain \code{"all"} mode).
#' @param max_tries     Positive integer. Maximum rejection-sampling attempts
#'   before falling back (default \code{100L}).
#'
#' @return Named list of sampled values.
#' @keywords internal
#' @noRd
.sample_from_space <- function(space,
                               gamma_ratio_min = NULL,
                               max_tries       = 100L) {

  space_names <- names(space)
  if (is.null(space_names) || length(space_names) == 0L)
    stop(".sample_from_space: `space` must be a *named* list.")
  space <- stats::setNames(
    lapply(space, function(spec) {
      if (!is.list(spec)) spec <- as.list(spec)
      if (!is.null(names(spec)) && "type" %in% names(spec)) return(spec)
      type <- spec[[1L]]
      switch(type,
             choice      = list(type = type, choices = spec[-1L]),
             randint     = list(type = type,
                                low  = as.integer(spec[[2L]]),
                                high = as.integer(spec[[3L]])),
             uniform     = list(type = type,
                                low  = as.numeric(spec[[2L]]),
                                high = as.numeric(spec[[3L]])),
             quniform    = list(type = type,
                                low  = as.numeric(spec[[2L]]),
                                high = as.numeric(spec[[3L]]),
                                q    = as.numeric(spec[[4L]])),
             loguniform  = list(type = type,
                                low  = as.numeric(spec[[2L]]),
                                high = as.numeric(spec[[3L]])),
             qloguniform = list(type = type,
                                low  = as.numeric(spec[[2L]]),
                                high = as.numeric(spec[[3L]]),
                                q    = as.numeric(spec[[4L]])),
             normal      = list(type = type,
                                mu    = as.numeric(spec[[2L]]),
                                sigma = as.numeric(spec[[3L]])),
             qnormal     = list(type = type,
                                mu    = as.numeric(spec[[2L]]),
                                sigma = as.numeric(spec[[3L]]),
                                q     = as.numeric(spec[[4L]])),
             lognormal   = list(type = type,
                                mu    = as.numeric(spec[[2L]]),
                                sigma = as.numeric(spec[[3L]])),
             qlognormal  = list(type = type,
                                mu    = as.numeric(spec[[2L]]),
                                sigma = as.numeric(spec[[3L]]),
                                q     = as.numeric(spec[[4L]])),
             stop(sprintf("Unknown search space type: '%s'", type))
      )
    }),
    space_names
  )

  # Sampler interne ---------------------------------------------------------
  .draw_once <- function() {
    params <- list()
    for (param_name in space_names) {
      spec <- space[[param_name]]
      params[[param_name]] <- switch(spec$type,
                                     choice      = sample(spec$choices, 1L)[[1L]],
                                     randint     = sample(seq.int(spec$low, spec$high), 1L),
                                     uniform     = stats::runif(1L, spec$low, spec$high),
                                     quniform    = {
                                       v <- stats::runif(1L, spec$low, spec$high)
                                       round(v / spec$q) * spec$q
                                     },
                                     loguniform  = exp(stats::runif(1L, log(spec$low), log(spec$high))),
                                     qloguniform = {
                                       v <- exp(stats::runif(1L, log(spec$low), log(spec$high)))
                                       round(v / spec$q) * spec$q
                                     },
                                     normal      = stats::rnorm(1L, spec$mu, spec$sigma),
                                     qnormal     = {
                                       v <- stats::rnorm(1L, spec$mu, spec$sigma)
                                       round(v / spec$q) * spec$q
                                     },
                                     lognormal   = stats::rlnorm(1L, spec$mu, spec$sigma),
                                     qlognormal  = {
                                       v <- stats::rlnorm(1L, spec$mu, spec$sigma)
                                       round(v / spec$q) * spec$q
                                     },
                                     stop(paste("Unknown search space type:", spec$type))
      )
    }
    if ("lambda_factor" %in% names(params) && "gamma" %in% names(params)) {
      params$lambda_ <- params$gamma * params$lambda_factor
    }

    params
  }
  if (is.null(gamma_ratio_min)) {
    return(.draw_once())
  }
  stopifnot(is.numeric(gamma_ratio_min), gamma_ratio_min > 0)

  params <- NULL
  for (attempt in seq_len(max_tries)) {
    candidate <- .draw_once()
    if (!("gamma" %in% names(candidate) && "lambda_" %in% names(candidate))) {
      return(candidate)
    }

    if (candidate$gamma > gamma_ratio_min * candidate$lambda_) {
      return(candidate)
    }
    params <- candidate
  }

  warning(sprintf(
    paste0(".sample_from_space: No satisfactory candidate found ",
           "Gamma > %.4g * lambda_ found in %d attempts.",
           "The last draw was returned as is"),
    gamma_ratio_min, max_tries
  ), call. = FALSE)

  params
}


# -----------------------------------------------------------------------------
# .custom_progress_bar  [private]
# -----------------------------------------------------------------------------

#' Minimal console progress bar
#'
#' @param total Positive integer. Total number of ticks.
#' @param width Integer. Bar width in characters (default \code{60}).
#' @return A list with one element \code{tick()}.
#' @keywords internal
#' @noRd
.custom_progress_bar <- function(total, width = 60L) {
  current <- 0L
  list(
    tick = function() {
      current  <<- current + 1L
      percent  <- current / total
      bars     <- round(width * percent)
      spaces   <- width - bars
      bar_text <- paste0(
        "[", paste(rep("=", bars),   collapse = ""),
        paste(rep(" ", spaces), collapse = ""),
        "] ", sprintf("%3.0f%%", percent * 100)
      )
      cat("\r", bar_text, sep = "", file = stderr())
      if (current == total) cat("\n", file = stderr())
    }
  )
}


# -----------------------------------------------------------------------------
# research_hyperOpt  [public]
# -----------------------------------------------------------------------------

#' Hyperparameter optimisation loop for dicepro
#'
#' Runs \code{hp_max_evals} NMF trials by sampling from \code{hp_space},
#' using either random or TPE sampling, and collects results.
#'
#' @param objective_opt Function passed verbatim to
#'   \code{objective_wrapper()}.
#' @param dataset       List with \code{$B}, \code{$W}, \code{$P}.
#' @param config        List. Configuration object; must contain
#'   \code{exp}, \code{hp_max_evals}, \code{hp_method}.
#'   Optionally contains \code{gamma_ratio_min} (positive numeric) to
#'   activate rejection sampling (mode \code{"all_gamma_dominant"}).
#' @param hp_space      Pre-built named list of parsed hyperparameter specs.
#'   If \code{NULL}, built from \code{config$hp_space}.
#' @param W_prime       Optional initial \eqn{W} matrix.
#' @param seed Integer. Random seed used for full pipeline reproducibility.
#'
#' @return A named list:
#' \describe{
#'   \item{trials}{data.frame of all successful trial results.}
#'   \item{W}{List of \eqn{W} matrices, one per successful trial.}
#'   \item{H}{List of \eqn{H} matrices, one per successful trial.}
#' }
#'
#' @export
research_hyperOpt <- function(objective_opt,
                              dataset,
                              config,
                              hp_space = NULL,
                              W_prime = NULL,
                              seed = NULL) {

  config <- .parse_config(config)

  seed <- seed %||% config$seed %||% 42L

  old_state <- if (exists(".Random.seed", .GlobalEnv)) .Random.seed else NULL

  set.seed(seed)

  on.exit({
    if (!is.null(old_state)) {
      .Random.seed <<- old_state
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  # Ratio minimal gamma/lambda_ pour le mode all_gamma_dominant (NULL = inactif)
  gamma_ratio_min <- config$gamma_ratio_min

  {
    search_space <- if (is.null(hp_space)) {
      stats::setNames(
        lapply(
          names(config$hp_space),
          function(arg) .parse_hyperopt_searchspace(arg, config$hp_space[[arg]])
        ),
        names(config$hp_space)
      )
    } else {
      hp_space
    }

    method  <- tolower(config$hp_method %||% "random")
    n_evals <- config$hp_max_evals

    results_list    <- vector("list", n_evals)
    W_list          <- vector("list", n_evals)
    out_deconv_list <- vector("list", n_evals)

    n_startup   <- max(10L, ceiling(sqrt(n_evals)))
    tpe_history <- list()

    pb <- .custom_progress_bar(total = n_evals, width = 60L)

    for (i in seq_len(n_evals)) {

      use_tpe <- method == "tpe" && length(tpe_history) >= n_startup

      params <- if (use_tpe) {
        .tpe_sample(search_space, tpe_history,
                    gamma_ratio_min = gamma_ratio_min)
      } else {
        .sample_from_space(search_space,
                           gamma_ratio_min = gamma_ratio_min)
      }

      res <- tryCatch(
        objective_wrapper(
          objective_opt = objective_opt,
          dataset       = dataset,
          config        = config,
          params        = params,
          W_prime       = W_prime
        ),
        error = function(e) {
          message(sprintf("Trial %d/%d failed: %s", i, n_evals, conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(res)) {
        results_list[[i]]    <- res$df
        W_list[[i]]          <- res$W
        out_deconv_list[[i]] <- res$H

        if (method == "tpe") {
          tpe_history <- c(tpe_history, list(list(
            params = params,
            loss   = res$df$loss
          )))
        }
      }

      pb$tick()
    }

    ok     <- !vapply(results_list, is.null, logical(1L))
    n_ok   <- sum(ok)
    n_fail <- n_evals - n_ok

    if (n_fail > 0L)
      message(sprintf("%d/%d trial(s) failed and were discarded.", n_fail, n_evals))

    if (n_ok == 0L) {
      warning("All ", n_evals, " trials failed -- returning empty results.")
      return(list(trials = data.frame(), W = list(), H = list()))
    }

    list(
      trials = do.call(rbind, results_list[ok]),
      W      = W_list[ok],
      H      = out_deconv_list[ok]
    )
  }
}
# -----------------------------------------------------------------------------
# .tpe_sample  [private]
# -----------------------------------------------------------------------------

#' TPE: sample one configuration guided by past trials
#'
#' @param search_space    Named list of parsed hyperparameter specs.
#' @param history         List of past trials (each a list with \code{params}
#'   and \code{loss}).
#' @param gamma           Numeric. Quantile threshold (default \code{0.25}).
#' @param n_candidates    Integer. Random candidates to score (default \code{24}).
#' @param gamma_ratio_min Numeric or \code{NULL}. Transmitted to
#'   \code{.sample_from_space()} for rejection sampling.
#'
#' @return Named list of sampled hyperparameter values.
#' @keywords internal
#' @noRd
.tpe_sample <- function(search_space, history,
                        gamma           = 0.25,
                        n_candidates    = 24L,
                        gamma_ratio_min = NULL) {

  losses    <- vapply(history, `[[`, numeric(1L), "loss")
  threshold <- stats::quantile(losses, probs = gamma, names = FALSE)

  good_idx <- which(losses <= threshold)
  bad_idx  <- which(losses >  threshold)

  if (length(good_idx) < 2L) good_idx <- order(losses)[seq_len(2L)]
  if (length(bad_idx)  < 2L) bad_idx  <- order(losses, decreasing = TRUE)[seq_len(2L)]

  candidates <- lapply(
    seq_len(n_candidates),
    function(.) .sample_from_space(search_space,
                                   gamma_ratio_min = gamma_ratio_min)
  )

  scores <- vapply(candidates, function(cand) {
    log_ratio <- 0
    for (pname in names(search_space)) {
      spec  <- search_space[[pname]]
      x_val <- cand[[pname]]
      if (is.null(x_val)) next

      good_vals <- vapply(history[good_idx],
                          function(h) h$params[[pname]] %||% NA_real_,
                          numeric(1L))
      bad_vals  <- vapply(history[bad_idx],
                          function(h) h$params[[pname]] %||% NA_real_,
                          numeric(1L))

      good_vals <- good_vals[!is.na(good_vals)]
      bad_vals  <- bad_vals[!is.na(bad_vals)]
      if (length(good_vals) < 2L || length(bad_vals) < 2L) next

      use_log <- spec$type %in%
        c("loguniform", "qloguniform", "lognormal", "qlognormal")

      if (use_log) {
        x_val     <- log(max(x_val,     .Machine$double.eps))
        good_vals <- log(pmax(good_vals, .Machine$double.eps))
        bad_vals  <- log(pmax(bad_vals,  .Machine$double.eps))
      }

      log_ratio <- log_ratio +
        .kde_log_density(x_val, bad_vals) -
        .kde_log_density(x_val, good_vals)
    }
    log_ratio
  }, numeric(1L))

  candidates[[which.max(scores)]]
}


# -----------------------------------------------------------------------------
# .kde_log_density  [private]
# -----------------------------------------------------------------------------

#' Gaussian KDE log-density at a single point
#'
#' Uses Silverman's rule of thumb for bandwidth selection.
#'
#' @param x   Numeric scalar. Query point.
#' @param obs Numeric vector. Observations.
#' @return Numeric scalar (log density).
#' @keywords internal
#' @noRd
.kde_log_density <- function(x, obs) {
  n      <- length(obs)
  sd_obs <- stats::sd(obs)

  bw <- if (sd_obs < .Machine$double.eps) 1e-3 else 1.06 * sd_obs * n^(-0.2)

  log_contribs <- stats::dnorm(x, mean = obs, sd = bw, log = TRUE)
  max_lc       <- max(log_contribs)
  log_sum      <- max_lc + log(sum(exp(log_contribs - max_lc)))

  log_sum - log(n)
}


# -----------------------------------------------------------------------------
# .parse_config  [private]
# -----------------------------------------------------------------------------

#' Validate and return a configuration list
#'
#' @param config Named list. Must contain \code{exp}, \code{hp_max_evals},
#'   and \code{hp_method}. Optionally contains \code{gamma_ratio_min}
#'   (positive numeric) used by the \code{"all_gamma_dominant"} mode.
#' @return The validated \code{config} list (unchanged).
#' @keywords internal
#' @noRd
.parse_config <- function(config) {
  required_args <- c("exp", "hp_max_evals", "hp_method")
  for (arg in required_args) {
    if (is.null(config[[arg]]))
      stop(paste("No", arg, "argument found in configuration file."),
           call. = FALSE)
  }

  valid_methods <- c("tpe", "random", "atpe", "anneal")
  if (!config$hp_method %in% valid_methods)
    stop(
      paste("Unknown hyperopt algorithm:", config$hp_method,
            "-- valid options:", paste(valid_methods, collapse = ", ")),
      call. = FALSE
    )

  # Validation of the optional ratio (all_gamma_dominant mode)
  if (!is.null(config$gamma_ratio_min)) {
    if (!is.numeric(config$gamma_ratio_min) ||
        length(config$gamma_ratio_min) != 1L ||
        config$gamma_ratio_min <= 0) {
      stop("'gamma_ratio_min' must be a single positive numeric.", call. = FALSE)
    }
  }

  config
}


# -----------------------------------------------------------------------------
# .parse_hyperopt_searchspace  [private]
# -----------------------------------------------------------------------------

#' Parse one hyperparameter search-space specification
#'
#' Converts an atomic vector \code{c(type, low, high, ...)} or an already
#' structured list into a fully named spec list consumed by
#' \code{.sample_from_space}.
#'
#' @param arg   Character scalar. Parameter name (used in error messages).
#' @param specs Atomic vector or list.
#' @return Named list with \code{type} and bounds/parameters.
#' @keywords internal
#' @noRd
.parse_hyperopt_searchspace <- function(arg, specs) {
  if (!is.list(specs)) specs <- as.list(specs)
  type <- specs[[1L]]
  .n   <- as.numeric

  switch(type,
         choice      = list(type = "choice",
                            choices = lapply(specs[-1L], .n)),
         randint     = list(type = "randint",
                            low  = as.integer(specs[[2L]]),
                            high = as.integer(specs[[3L]])),
         uniform     = list(type = "uniform",
                            low  = .n(specs[[2L]]),
                            high = .n(specs[[3L]])),
         quniform    = list(type = "quniform",
                            low  = .n(specs[[2L]]),
                            high = .n(specs[[3L]]),
                            q    = .n(specs[[4L]])),
         loguniform  = list(type = "loguniform",
                            low  = .n(specs[[2L]]),
                            high = .n(specs[[3L]])),
         qloguniform = list(type = "qloguniform",
                            low  = .n(specs[[2L]]),
                            high = .n(specs[[3L]]),
                            q    = .n(specs[[4L]])),
         normal      = list(type = "normal",
                            mu    = .n(specs[[2L]]),
                            sigma = .n(specs[[3L]])),
         qnormal     = list(type = "qnormal",
                            mu    = .n(specs[[2L]]),
                            sigma = .n(specs[[3L]]),
                            q     = .n(specs[[4L]])),
         lognormal   = list(type = "lognormal",
                            mu    = .n(specs[[2L]]),
                            sigma = .n(specs[[3L]])),
         qlognormal  = list(type = "qlognormal",
                            mu    = .n(specs[[2L]]),
                            sigma = .n(specs[[3L]]),
                            q     = .n(specs[[4L]])),
         stop(sprintf(
           "Unknown search space type '%s' for parameter '%s'.", type, arg
         ), call. = FALSE)
  )
}


# -----------------------------------------------------------------------------
# .get_report_path  [private]
# -----------------------------------------------------------------------------

#' Ensure a results subdirectory exists and return its path
#'
#' @param exp_dir Character scalar. Experiment output directory.
#' @return Character scalar. Path to \code{exp_dir/results/}.
#' @keywords internal
#' @noRd
.get_report_path <- function(exp_dir) {
  report_path <- file.path(exp_dir, "results")
  if (!dir.exists(report_path))
    dir.create(report_path, recursive = TRUE, showWarnings = FALSE)
  report_path
}
