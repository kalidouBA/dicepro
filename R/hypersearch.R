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

# ---------------------------------------------------------------------------
# Internal helper: clean an H or W matrix returned by nmf_lbfgsb.
# ---------------------------------------------------------------------------
.clean_nmf_matrix <- function(mat) {
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (!is.matrix(mat)) return(mat)

  char_cols <- which(apply(mat, 2, function(x) !suppressWarnings(all(!is.na(as.numeric(x))))))

  if (length(char_cols) > 0L) {
    rn      <- mat[, char_cols[1L]]
    mat     <- mat[, -char_cols, drop = FALSE]
    rownames(mat) <- rn
  }

  rn <- rownames(mat)
  cn <- colnames(mat)
  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat),
                dimnames = list(rn, cn))
  mat
}

#' NMF LBFGSB Hyperparameter Optimization
#'
#' @param dataset List containing matrices `B`, `W`, and `P`.
#' @param W_prime Optional numeric matrix. Initial W'.
#' @param p_prime Optional numeric matrix. Initial P'.
#' @param lambda_ Numeric. Regularization parameter lambda.
#' @param gamma_par Numeric. Regularization parameter gamma.
#' @param path2save Character. Path to save results.
#' @return List. Output from `nmf_lbfgsb`.
#' @import dicepro
#' @export
nmf_lbfgsb_hyperOpt <- function(dataset, W_prime = NULL, p_prime = NULL,
                                lambda_ = 10, gamma_par = 100, path2save = "") {
  B    <- as.data.frame(dataset$B)
  W_cb <- as.data.frame(dataset$W)
  P_cb <- as.data.frame(dataset$P)

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
  }

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

#' Objective function for hyperparameter optimization
#'
#' @param dataset List. Dataset containing matrices `B`, `W`, and `P`.
#' @param config List. Configuration for the experiment.
#' @param lambda_ Numeric. Regularization parameter lambda.
#' @param gamma_factor Numeric or NULL. If provided: gamma = lambda_ * gamma_factor.
#' @param gamma Numeric. Regularization parameter gamma (used directly when gamma_factor is NULL).
#' @param p_prime Numeric matrix (nrow = n_samples, ncol = n_unknown_ct).
#' @param W_prime Numeric. Optional matrix for initialization.
#' @return List with loss, constraint, status, current_params, W, H. NULL if invalid.
#' @export
objective_opt <- function(dataset, config = list(), lambda_ = NULL,
                          gamma_factor = NULL, gamma = NULL,
                          p_prime = NULL, W_prime = 0) {

  exp_dir <- ifelse(!is.null(config$exp), config$exp, ".")

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

  result_dict <- lapply(names(result), function(nm) {
    x <- result[[nm]]
    if (nm %in% c("H", "W", "p_prime_estm")) {
      if (is.data.frame(x)) as.matrix(x) else x
    } else if (is.numeric(x) && (is.null(dim(x)) || length(dim(x)) == 0)) {
      x
    } else if (is.matrix(x) || is.data.frame(x)) {
      num_cols <- sapply(x, is.numeric)
      as.numeric(as.matrix(x[, num_cols, drop = FALSE]))
    } else {
      x
    }
  })
  names(result_dict) <- names(result)

  .sc <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    as.numeric(x[[1L]])
  }

  if (contains_nan_or_inf(.sc(result_dict$objectiveValue)) ||
      contains_nan_or_inf(.sc(result_dict$constraint)) ||
      any(is.na(c(.sc(result_dict$frob_H), .sc(result_dict$var_H),
                  .sc(result_dict$frob_W), .sc(result_dict$var_W))))) {
    return(NULL)
  }

  return(list(
    loss       = .sc(result_dict$objectiveValue),
    constraint = .sc(result_dict$constraint),
    status     = "OK",
    current_params = list(
      gamma          = gamma,
      lambda_        = lambda_,
      p_prime        = p_prime[1, 1],
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
  ))
}

#' Objective wrapper for hyperparameter optimization
#'
#' @param objective_opt Function. The objective function to optimize.
#' @param dataset List. Dataset used for optimization.
#' @param config List. Parsed configuration object.
#' @param params List. Must contain `lambda_`, `gamma`, `p_prime`; optionally `gamma_factor`.
#' @param W_prime Matrix or scalar. Optional initial W matrix for NMF.
#' @return list(df, W, H) or NULL if trial failed.
#' @importFrom utils capture.output
#' @keywords internal
objective_wrapper <- function(objective_opt, dataset, config, params, W_prime = NULL) {

  tryCatch(
    {
      start_time <- Sys.time()

      n_samples    <- ncol(dataset$B)
      n_unknown_ct <- if (!is.null(W_prime) && is.matrix(W_prime) && ncol(W_prime) > 0L)
        ncol(W_prime) else 1L

      if (!is.matrix(params$p_prime)) {
        params$p_prime <- matrix(params$p_prime, nrow = n_samples, ncol = n_unknown_ct)
      }

      gamma_factor_arg <- params$gamma_factor

      suppressMessages(
        capture.output(
          invisible(returned_dict <- objective_opt(
            dataset      = dataset,
            config       = config,
            lambda_      = params$lambda_,
            gamma_factor = gamma_factor_arg,
            gamma        = params$gamma,
            p_prime      = params$p_prime,
            W_prime      = W_prime
          )),
          file = "/dev/null"
        )
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

      list(df = as.data.frame(row, check.names = FALSE),
           W  = returned_dict$W,
           H  = returned_dict$H)

    },
    error = function(e) {
      message("objective_wrapper error: ", conditionMessage(e))
      return(NULL)
    })
}


#' Sample from hyperparameter space
#'
#' @param space Named list. Each element is a named list with `type` and bounds.
#' @return List of sampled parameter values.
#' @noRd
.sample_from_space <- function(space) {
  space_names <- names(space)
  if (is.null(space_names) || length(space_names) == 0L)
    stop(".sample_from_space: `space` must be a *named* list.")

  space <- stats::setNames(
    lapply(space, function(spec) {
      if (!is.list(spec)) spec <- as.list(spec)
      if (!is.null(names(spec)) && "type" %in% names(spec)) return(spec)
      type <- spec[[1]]
      switch(type,
             "choice"      = list(type = type, choices = spec[-1]),
             "randint"     = list(type = type,
                                  low  = as.integer(spec[[2]]),
                                  high = as.integer(spec[[3]])),
             "uniform"     = list(type = type,
                                  low  = as.numeric(spec[[2]]),
                                  high = as.numeric(spec[[3]])),
             "quniform"    = list(type = type,
                                  low  = as.numeric(spec[[2]]),
                                  high = as.numeric(spec[[3]]),
                                  q    = as.numeric(spec[[4]])),
             "loguniform"  = list(type = type,
                                  low  = as.numeric(spec[[2]]),
                                  high = as.numeric(spec[[3]])),
             "qloguniform" = list(type = type,
                                  low  = as.numeric(spec[[2]]),
                                  high = as.numeric(spec[[3]]),
                                  q    = as.numeric(spec[[4]])),
             "normal"      = list(type = type,
                                  mu    = as.numeric(spec[[2]]),
                                  sigma = as.numeric(spec[[3]])),
             "qnormal"     = list(type = type,
                                  mu    = as.numeric(spec[[2]]),
                                  sigma = as.numeric(spec[[3]]),
                                  q     = as.numeric(spec[[4]])),
             "lognormal"   = list(type = type,
                                  mu    = as.numeric(spec[[2]]),
                                  sigma = as.numeric(spec[[3]])),
             "qlognormal"  = list(type = type,
                                  mu    = as.numeric(spec[[2]]),
                                  sigma = as.numeric(spec[[3]]),
                                  q     = as.numeric(spec[[4]])),
             stop(sprintf("Unknown search space type: '%s'", type))
      )
    }),
    space_names
  )

  params <- list()

  for (param_name in space_names) {
    spec <- space[[param_name]]

    params[[param_name]] <- switch(spec$type,
                                   "choice"      = sample(spec$choices, 1L)[[1L]],
                                   "randint"     = sample(seq.int(spec$low, spec$high), 1L),
                                   "uniform"     = runif(1L, spec$low, spec$high),
                                   "quniform"    = {
                                     v <- runif(1L, spec$low, spec$high)
                                     round(v / spec$q) * spec$q
                                   },
                                   "loguniform"  = exp(runif(1L, log(spec$low), log(spec$high))),
                                   "qloguniform" = {
                                     v <- exp(runif(1L, log(spec$low), log(spec$high)))
                                     round(v / spec$q) * spec$q
                                   },
                                   "normal"      = rnorm(1L, spec$mu, spec$sigma),
                                   "qnormal"     = {
                                     v <- rnorm(1L, spec$mu, spec$sigma)
                                     round(v / spec$q) * spec$q
                                   },
                                   "lognormal"   = rlnorm(1L, spec$mu, spec$sigma),
                                   "qlognormal"  = {
                                     v <- rlnorm(1L, spec$mu, spec$sigma)
                                     round(v / spec$q) * spec$q
                                   },
                                   stop(paste("Unknown search space type:", spec$type))
    )
  }

  if ("lambda_factor" %in% names(params) && "gamma" %in% names(params)) {
    params$lambda_ <- params$gamma * params$lambda_factor
  }

  return(params)
}


#' Custom progress bar
#' @noRd
.custom_progress_bar <- function(total, format = "Progress [:bar] :percent", width = 60) {
  current <- 0
  list(
    tick = function() {
      current <<- current + 1
      percent  <- current / total
      bars     <- round(width * percent)
      spaces   <- width - bars
      bar_text <- paste0("[", paste(rep("=", bars), collapse = ""),
                         paste(rep(" ", spaces), collapse = ""),
                         "] ", sprintf("%3.0f%%", percent * 100))
      message("\r", bar_text, appendLF = FALSE)
      if (current == total) message()
    }
  )
}


#' Hyperparameter optimisation for dicepro
#'
#' @param objective_opt  Function passed verbatim to \code{objective_wrapper()}.
#' @param dataset        List with \code{$B}, \code{$W}, \code{$P}.
#' @param config         List. Configuration object (replaces config_path).
#' @param hp_space       Pre-built space list; if NULL, built from config$hp_space.
#' @param W_prime        Optional initial W matrix.
#' @return list(trials, W, H).
#' @export
research_hyperOpt <- function(objective_opt, dataset,
                              config,                # <-- plus config_path
                              hp_space = NULL,
                              W_prime = NULL) {
  config <- .parse_config(config)
  set.seed(config$seed %||% 42L)

  search_space <- if (is.null(hp_space)) {
    stats::setNames(
      lapply(names(config$hp_space), function(arg)
        .parse_hyperopt_searchspace(arg, config$hp_space[[arg]])
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
      .tpe_sample(search_space, tpe_history)
    } else {
      .sample_from_space(search_space)
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


#' TPE: sample one set of hyperparameters guided by past trials
#'
#' @param search_space Named list of parsed hyperparameter specs.
#' @param history List of past trials.
#' @param gamma Quantile threshold (default 0.25).
#' @param n_candidates Number of random candidates to score (default 24).
#' @return Named list of sampled hyperparameter values.
#' @noRd
.tpe_sample <- function(search_space, history,
                        gamma        = 0.25,
                        n_candidates = 24L) {

  losses <- vapply(history, `[[`, numeric(1L), "loss")

  threshold <- stats::quantile(losses, probs = gamma, names = FALSE)

  good_idx <- which(losses <= threshold)
  bad_idx  <- which(losses >  threshold)

  if (length(good_idx) < 2L) good_idx <- order(losses)[seq_len(2L)]
  if (length(bad_idx)  < 2L) bad_idx  <- order(losses, decreasing = TRUE)[seq_len(2L)]

  candidates <- lapply(seq_len(n_candidates), function(.) .sample_from_space(search_space))

  scores <- vapply(candidates, function(cand) {
    log_ratio <- 0
    for (pname in names(search_space)) {
      spec  <- search_space[[pname]]
      x_val <- cand[[pname]]
      if (is.null(x_val)) next

      good_vals <- vapply(history[good_idx], function(h) h$params[[pname]] %||% NA_real_,
                          numeric(1L))
      bad_vals  <- vapply(history[bad_idx],  function(h) h$params[[pname]] %||% NA_real_,
                          numeric(1L))

      good_vals <- good_vals[!is.na(good_vals)]
      bad_vals  <- bad_vals[!is.na(bad_vals)]
      if (length(good_vals) < 2L || length(bad_vals) < 2L) next

      use_log <- spec$type %in% c("loguniform", "qloguniform", "lognormal", "qlognormal")
      if (use_log) {
        x_val     <- log(max(x_val,     .Machine$double.eps))
        good_vals <- log(pmax(good_vals, .Machine$double.eps))
        bad_vals  <- log(pmax(bad_vals,  .Machine$double.eps))
      }

      log_p_good <- .kde_log_density(x_val, good_vals)
      log_p_bad  <- .kde_log_density(x_val, bad_vals)

      log_ratio <- log_ratio + (log_p_bad - log_p_good)
    }
    log_ratio
  }, numeric(1L))

  best_idx <- which.max(scores)
  candidates[[best_idx]]
}


#' Gaussian KDE log-density at a single point
#'
#' @param x   Scalar query point.
#' @param obs Numeric vector of observations.
#' @return Numeric scalar (log density).
#' @noRd
.kde_log_density <- function(x, obs) {
  n      <- length(obs)
  sd_obs <- stats::sd(obs)

  if (sd_obs < .Machine$double.eps) {
    bw <- 1e-3
  } else {
    bw <- 1.06 * sd_obs * n^(-0.2)
  }

  log_contribs <- stats::dnorm(x, mean = obs, sd = bw, log = TRUE)

  max_lc  <- max(log_contribs)
  log_sum <- max_lc + log(sum(exp(log_contribs - max_lc)))

  log_sum - log(n)
}


#' Parse configuration list
#' @noRd
.parse_config <- function(config) {
  required_args <- c("exp", "hp_max_evals", "hp_method")
  for (arg in required_args) {
    if (is.null(config[[arg]]))
      stop(paste("No", arg, "argument found in configuration file."))
  }
  valid_methods <- c("tpe", "random", "atpe", "anneal")
  if (!config$hp_method %in% valid_methods)
    stop(paste("Unknown hyperopt algorithm:", config$hp_method,
               "-- valid options:", paste(valid_methods, collapse = ", ")))
  return(config)
}

#' Parse hyperopt search space specification
#' @noRd
.parse_hyperopt_searchspace <- function(arg, specs) {
  if (!is.list(specs)) specs <- as.list(specs)
  type <- specs[[1]]
  .n   <- function(x) as.numeric(x)

  switch(type,
         "choice"      = list(type = "choice",      choices = lapply(specs[-1], .n)),
         "randint"     = list(type = "randint",      low = as.integer(specs[[2]]),
                              high = as.integer(specs[[3]])),
         "uniform"     = list(type = "uniform",      low = .n(specs[[2]]), high = .n(specs[[3]])),
         "quniform"    = list(type = "quniform",     low = .n(specs[[2]]), high = .n(specs[[3]]),
                              q = .n(specs[[4]])),
         "loguniform"  = list(type = "loguniform",   low = .n(specs[[2]]), high = .n(specs[[3]])),
         "qloguniform" = list(type = "qloguniform",  low = .n(specs[[2]]), high = .n(specs[[3]]),
                              q = .n(specs[[4]])),
         "normal"      = list(type = "normal",       mu = .n(specs[[2]]), sigma = .n(specs[[3]])),
         "qnormal"     = list(type = "qnormal",      mu = .n(specs[[2]]), sigma = .n(specs[[3]]),
                              q = .n(specs[[4]])),
         "lognormal"   = list(type = "lognormal",    mu = .n(specs[[2]]), sigma = .n(specs[[3]])),
         "qlognormal"  = list(type = "qlognormal",   mu = .n(specs[[2]]), sigma = .n(specs[[3]]),
                              q = .n(specs[[4]])),
         stop(sprintf("Unknown search space type '%s' for parameter '%s'.", type, arg))
  )
}


#' Get report path
#' @noRd
.get_report_path <- function(exp_dir) {
  report_path <- file.path(exp_dir, "results")
  if (!dir.exists(report_path)) dir.create(report_path, recursive = TRUE)
  return(report_path)
}
