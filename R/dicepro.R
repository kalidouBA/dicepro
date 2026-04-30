# =============================================================================
# dicepro.R -- Main entry point
#
# Public API  : dicepro()
# Private fns : .validate_inputs(), .normalize_zscore_per_gene(),
#               .prepare_data(), .get_deconvolution(),
#               .run_hyperopt(), .save_outputs()
# =============================================================================


# -----------------------------------------------------------------------------
# .validate_inputs  [private]
# -----------------------------------------------------------------------------

#' Validate and normalize dicepro user parameters
#'
#' @param normalize             Logical.
#' @param algo_select           Character. hyper-parameter search strategy.
#' @param hspaceTechniqueChoose Character. Search-space strategy.
#' @param gamma_ratio_min       Numeric or \code{NULL}.
#' @param constraint_threshold Positive numeric.
#' @param cibersort_perm        Non-negative integer.
#' @param cibersort_QN          Logical.
#'
#' @return Named list with validated fields.
#'
#' @keywords internal
#' @noRd
.validate_inputs <- function(normalize,
                             algo_select,
                             hspaceTechniqueChoose,
                             gamma_ratio_min = NULL,
                             constraint_threshold = 0.1,
                             cibersort_perm  = 0L,
                             cibersort_QN    = TRUE) {

  if (!is.logical(normalize) || length(normalize) != 1L)
    stop("'normalize' must be TRUE or FALSE.", call. = FALSE)

  algo_select <- match.arg(
    tolower(algo_select),
    c("random", "tpe", "atpe", "anneal")
  )

  hspaceTechniqueChoose <- match.arg(
    hspaceTechniqueChoose,
    c("all", "restrictionEspace", "gamma_dominant")
  )

  if (hspaceTechniqueChoose == "gamma_dominant") {
    if (is.null(gamma_ratio_min)) gamma_ratio_min <- 10
    if (!is.numeric(gamma_ratio_min) ||
        length(gamma_ratio_min) != 1L ||
        gamma_ratio_min <= 0)
      stop("'gamma_ratio_min' must be a single positive numeric.", call. = FALSE)
  }

  if (!is.numeric(constraint_threshold) ||
      length(constraint_threshold) != 1L ||
      is.na(constraint_threshold) ||
      !is.finite(constraint_threshold) ||
      constraint_threshold <= 0) {
    stop(
      "'constraint_threshold' must be a single positive finite numeric value.\n",
      "Guidelines: 0.05 = very strict (risk of no solution), ",
      "0.1 = good trade-off, 0.2 -- 0.3 = more permissive.",
      call. = FALSE
    )
  }

  # ---- CIBERSORT-specific validation ---------------------------------------
  if (!is.numeric(cibersort_perm) ||
      length(cibersort_perm) != 1L ||
      cibersort_perm < 0)
    stop("'cibersort_perm' must be a single non-negative integer.", call. = FALSE)
  cibersort_perm <- as.integer(cibersort_perm)

  if (!is.logical(cibersort_QN) || length(cibersort_QN) != 1L)
    stop("'cibersort_QN' must be TRUE or FALSE.", call. = FALSE)

  list(
    normalize             = normalize,
    algo_select           = algo_select,
    hspaceTechniqueChoose = hspaceTechniqueChoose,
    gamma_ratio_min       = gamma_ratio_min,
    constraint_threshold  = constraint_threshold,  # <-- IMPORTANT
    cibersort_perm        = cibersort_perm,
    cibersort_QN          = cibersort_QN
  )
}

# -----------------------------------------------------------------------------
# .normalize_zscore_per_gene  [private]
# -----------------------------------------------------------------------------

#' Z-score normalization per gene
#'
#' Center and scales each row (gene) of a numeric matrix. Genes with
#' \code{SD = 0} or \code{NA} are silently removed.
#'
#' @param mat Numeric matrix (genes * samples).
#' @return Numeric matrix with mean = 0 and SD = 1 per row.
#'
#' @keywords internal
#' @noRd
.normalize_zscore_per_gene <- function(mat) {

  mat     <- as.matrix(mat)
  gene_sd <- apply(mat, 1L, stats::sd, na.rm = TRUE)
  keep    <- !is.na(gene_sd) & gene_sd > 0
  mat     <- mat[keep, , drop = FALSE]

  t(apply(mat, 1L, function(x) {
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  }))
}


# -----------------------------------------------------------------------------
# .prepare_data  [private]
# -----------------------------------------------------------------------------

#' normalize and intersect bulk and reference matrices
#'
#' @param reference Numeric matrix (genes * cell types).
#' @param bulk      Numeric matrix (genes * samples).
#' @param normalize Logical.
#'
#' @return Named list: \code{reference}, \code{bulk} -- both filtered to
#'   intersecting genes.
#'
#' @keywords internal
#' @noRd
.prepare_data <- function(reference, bulk, normalize) {

  if (normalize) {
    reference <- .normalize_zscore_per_gene(reference)
    bulk      <- .normalize_zscore_per_gene(bulk)
  } else {
    reference <- as.matrix(reference)
    bulk      <- as.matrix(bulk)
  }

  geneIntersect <- intersect(rownames(reference), rownames(bulk))

  if (length(geneIntersect) == 0L)
    stop("No common genes between 'reference' and 'bulk'.", call. = FALSE)

  reference <- reference[geneIntersect, , drop = FALSE]
  bulk      <- bulk[geneIntersect, , drop = FALSE]
  rownames(reference) <- rownames(bulk) <- geneIntersect

  message(sprintf("Gene intersection: %d genes retained.", length(geneIntersect)))

  list(reference = reference, bulk = bulk)
}


# -----------------------------------------------------------------------------
# .get_deconvolution  [private]
# -----------------------------------------------------------------------------

#' Run or validate supervised deconvolution
#'
#' Calls \code{\link{running_method}} when \code{out_Decon} is \code{NULL};
#' otherwise validates and transposes the user-supplied matrix.
#'
#' @param bulk,reference         Numeric matrices (genes * sample * N_cell_type).
#' @param methodDeconv           Character. Deconvolution method.
#' @param cibersortx_email,cibersortx_token CIBERSORTx credentials.
#' @param out_Decon              Optional pre-computed matrix
#'   (samples * cell types) or \code{NULL}.
#' @param cibersort_perm         Non-negative integer. Permutations for
#'   built-in CIBERSORT p-values.
#' @param cibersort_QN           Logical. Quantile-normalize mixture in
#'   built-in CIBERSORT.
#'
#' @return Numeric matrix (cell types * samples).
#'
#' @keywords internal
#' @noRd
.get_deconvolution <- function(bulk,
                               reference,
                               methodDeconv,
                               cibersortx_email,
                               cibersortx_token,
                               out_Decon,
                               cibersort_perm = 0L,
                               cibersort_QN   = TRUE) {

  valid_methods <- c("CSx", "CS", "DCQ", "FARDEEP")

  if (is.null(out_Decon)) {

    methodDeconv <- tryCatch(
      match.arg(methodDeconv, valid_methods),
      error = function(e) stop(sprintf(
        "Invalid 'methodDeconv': '%s'. Must be one of: %s.",
        methodDeconv, paste(valid_methods, collapse = ", ")
      ), call. = FALSE)
    )

    if (methodDeconv == "CSx" &&
        (is.null(cibersortx_email) || is.null(cibersortx_token)))
      stop("CIBERSORTx credentials required for methodDeconv = 'CSx'.",
           call. = FALSE)

    out_Decon <- running_method(
      bulk             = bulk,
      reference        = reference,
      methodDeconv     = methodDeconv,
      cibersortx_email = cibersortx_email,
      cibersortx_token = cibersortx_token,
      cibersort_perm   = cibersort_perm,
      cibersort_QN     = cibersort_QN
    )

  } else {

    out_Decon <- as.matrix(out_Decon)

    if (ncol(out_Decon) != ncol(reference))
      stop("Mismatch: number of cell types vs reference columns.", call. = FALSE)

    if (nrow(out_Decon) != ncol(bulk))
      stop("Mismatch: number of samples vs bulk columns.", call. = FALSE)
  }

  t(out_Decon)
}


# -----------------------------------------------------------------------------
# .run_hyperopt  [private]
# -----------------------------------------------------------------------------

#' Execute the full hyper-parameter optimization pipeline
#'
#' @param dataset               Named list: \code{B}, \code{W}, \code{P}.
#' @param W_prime               Numeric matrix or scalar.
#' @param bulkName,refName      Character scalars.
#' @param hp_max_evals          Positive integer.
#' @param algo_select           Character scalar.
#' @param output_base_dir       Character scalar.
#' @param hspaceTechniqueChoose Character scalar.
#' @param output_dir            Character scalar.
#' @param gamma_ratio_min       Numeric or \code{NULL}.
#' @param seed                  Integer. Random seed.
#'
#' @return Output of \code{\link{best_hyperParams}}, or \code{NULL}.
#'
#' @keywords internal
#' @noRd
.run_hyperopt <- function(dataset,
                          W_prime,
                          bulkName,
                          refName,
                          hp_max_evals,
                          algo_select,
                          output_base_dir,
                          hspaceTechniqueChoose,
                          output_dir,
                          gamma_ratio_min,
                          constraint_threshold,
                          seed) {

  res <- run_experiment(
    dataset               = dataset,
    W_prime               = W_prime,
    bulkName              = bulkName,
    refName               = refName,
    hp_max_evals          = hp_max_evals,
    algo_select           = algo_select,
    output_base_dir       = output_base_dir,
    hspaceTechniqueChoose = hspaceTechniqueChoose,
    gamma_ratio_min       = gamma_ratio_min,
    seed                  = seed
  )

  best_hyperParams(
    trials_df = res$trials,
    W         = res$W,
    H         = res$H,
    savePaths = output_dir,
    constraint_threshold = constraint_threshold
  )
}


# -----------------------------------------------------------------------------
# .save_outputs  [private]
# -----------------------------------------------------------------------------

#' Save hyper-parameter optimization diagnostic plots
#'
#' Creates \code{output_dir/report/}, saves two PDFs (hyperopt report and
#' Pareto frontier), and attaches the hyperopt plot to \code{out}.
#'
#' @param out        Named list returned by \code{.run_hyperopt}.
#' @param output_dir Character scalar. Root results directory.
#' @param hp_params  Character vector of hyper-parameter names to plot.
#'
#' @return \code{out} with an additional \code{plot_hyperopt} element.
#'
#' @keywords internal
#' @noRd
.save_outputs <- function(out, output_dir, hp_params) {

  report_dir <- file.path(output_dir, "report")
  if (!dir.exists(report_dir))
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  out$plot_hyperopt <- plot_hyperopt(
    x      = structure(out, class = "dicepro"),
    params = hp_params
  )

  ggplot2::ggsave(
    filename = file.path(report_dir, "hyperopt_report.pdf"),
    plot     = out$plot_hyperopt,
    width    = 4L * length(hp_params),
    height   = 4L * length(hp_params) + 2L,
    device   = grDevices::cairo_pdf
  )

  ggplot2::ggsave(
    filename = file.path(report_dir, "pareto_frontier.pdf"),
    plot     = out$plot,
    width    = 8L,
    height   = 6L,
    device   = grDevices::cairo_pdf
  )

  message(sprintf("Results saved to: %s", report_dir))
  out
}


# -----------------------------------------------------------------------------
# dicepro  [public]
# -----------------------------------------------------------------------------

#' Semi-supervised bulk RNA-seq deconvolution with hyper-parameter optimization
#'
#' @description
#' Combines supervised estimation of known cell types with unsupervised
#' discovery of latent components, using automatic Pareto-frontier-based
#' hyper-parameter optimization.
#'
#' @details
#' When \code{out_Decon} is provided, the supervised step is skipped.
#' Gene matrices are optionally z-score normalized per gene, and only
#' intersecting genes are retained before optimization.
#'
#' The built-in \code{"CIBERSORT"} method closely follows Newman et al.
#' (2015): nu-SVR with three candidate nu values (0.25, 0.50, 0.75),
#' optional quantile normalization, and optional permutation-based p-values.
#' It requires packages \pkg{e1071}, \pkg{parallel}, and
#' \pkg{preprocessCore}. No external account or Docker installation is
#' needed, unlike \code{"CSx"}.
#'
#' @param reference             Numeric matrix (genes * cell types).
#' @param bulk                  Numeric matrix (genes * samples).
#' @param methodDeconv          Character. One of \code{"CSx"},
#'   \code{"CS"}, \code{"DCQ"}, \code{"FARDEEP"}.
#' @param cibersortx_email      Character. CIBERSORTx email (required for
#'   \code{"CSx"}).
#' @param cibersortx_token      Character. CIBERSORTx token (required for
#'   \code{"CSx"}).
#' @param cibersort_perm        Non-negative integer. Number of permutations
#'   for p-value estimation in the built-in \code{"CIBERSORT"} method.
#'   Set to \code{0L} (default) to skip p-value computation. Ignored for
#'   all other methods.
#' @param cibersort_QN          Logical. Apply quantile normalization to the
#'   bulk mixture in the built-in \code{"CIBERSORT"} method (default
#'   \code{TRUE}). Ignored for all other methods.
#' @param W_prime               Initial unknown-signature matrix or \code{0}.
#' @param bulkName              Character scalar. Label for the bulk data-set.
#' @param refName               Character scalar. Label for the reference.
#' @param hp_max_evals          Positive integer. Number of hyper-parameter
#'   trials.
#' @param N_unknownCT           Positive integer. Number of unknown cell
#'   types.
#' @param algo_select           Character. One of \code{"random"},
#'   \code{"tpe"}, \code{"atpe"}, \code{"anneal"}.
#' @param output_path           Character scalar. Root output directory
#'   (\code{getwd()} by default).
#' @param hspaceTechniqueChoose Character. One of \code{"gamma_dominant"}
#'  or \code{"all"}, or \code{"restrictionEspace"}.
#'   \strong{The most efficient and strongly recommended strategy is
#'   \code{"gamma_dominant"}, as it provides faster convergence and
#'   more precise results.}
#'   \describe{
#'     \item{\code{"gamma_dominant"}}{\strong{Recommended.}
#'       Same unconstrained grid as \code{"all"}, but candidates where
#'       \code{gamma <= gamma_ratio_min * lambda_} are rejected, ensuring
#'       \code{gamma >> lambda_}. This reduces the search space, leading
#'       to faster optimization and more stable solutions.}
#'
#'     \item{\code{"all"}}{Full independent log-uniform grid for
#'       \code{lambda_}, \code{gamma}, \code{p_prime}.}
#'
#'     \item{\code{"restrictionEspace"}}{Restricted space where
#'       \code{lambda_ = gamma * lambda_factor},
#'       \code{lambda_factor} \eqn{\in [2, 100]}.}
#'   }
#' @param gamma_ratio_min       Positive numeric. Minimum ratio
#'   \eqn{\gamma / \lambda} enforced when
#'   \code{hspaceTechniqueChoose = "gamma_dominant"} (default \code{10}).
#'   Ignored for other strategies.
#' @param constraint_threshold Positive numeric. Maximum allowed deviation
#'   from the constraint (|1 - constraint|) when filtering hyper-parameter
#'   trials before Pareto selection (default \code{0.1}).
#'   \describe{
#'     \item{\code{0.05}}{Very strict (risk of no valid solution).}
#'     \item{\code{0.1}}{Good trade-off between feasibility and robustness.}
#'     \item{\code{0.2 -- 0.3}}{More permissive, useful for broader exploration.}
#'   }
#' @param out_Decon             Optional pre-computed deconvolution matrix
#'   (samples * cell types). When provided, the supervised deconvolution
#'   step is skipped entirely.
#' @param normalize             Logical. Apply z-score normalization per
#'   gene (default \code{TRUE}).
#' @param seed                  Integer. Global random seed for
#'   reproducibility (default \code{42L}).
#'
#' @return An object of class \code{"dicepro"} (a named list) containing:
#' \itemize{
#'   \item \code{hyper-parameters} -- selected \eqn{\lambda} and \eqn{\gamma}
#'   \item \code{metrics}         -- loss and constraint of the best trial
#'   \item \code{trials}          -- all evaluated configurations
#'   \item \code{W}               -- estimated unknown signature matrix
#'   \item \code{H}               -- estimated proportion matrix
#'   \item \code{plot}            -- Pareto frontier ggplot2 figure
#'   \item \code{plot_hyperopt}   -- hyper-parameter space ggplot2 figure
#' }
#' Returns \code{invisible(NULL)} with a warning when no valid configuration
#' is found.
#'
#' @seealso \code{\link{running_method}}, \code{\link{run_experiment}},
#'   \code{\link{best_hyperParams}}
#'
#' @references
#' Newman AM et al. (2015). Robust enumeration of cell subsets from tissue
#' expression profiles. \emph{Nature Methods}, 12(5), 453-457.
#' \doi{10.1038/nmeth.3337}
#'
#' @examples
#' \dontrun{
#'  # gamma >> lambda (gamma_dominant strategy)
#'     # FARDEEP deconvolution
#' res <- dicepro(
#'   reference    = BlueCode,
#'   bulk         = CellMixtures,
#'   methodDeconv = "FARDEEP",
#'   hp_max_evals = 50L
#' )
#'
#'     # Built-in CIBERSORT (no Docker, no account)
#' res <- dicepro(
#'   reference      = BlueCode,
#'   bulk           = CellMixtures,
#'   methodDeconv   = "CS",
#'   cibersort_perm = 100L,
#'   cibersort_QN   = TRUE,
#'   hp_max_evals   = 50L
#' )
#' }
#'
#' @export
dicepro <- function(reference,
                    bulk,
                    methodDeconv          = "CSx",
                    cibersortx_email      = NULL,
                    cibersortx_token      = NULL,
                    cibersort_perm        = 0L,
                    cibersort_QN          = TRUE,
                    W_prime               = 0,
                    bulkName              = "Bulk",
                    refName               = "Reference",
                    hp_max_evals          = 100L,
                    N_unknownCT           = 1L,
                    algo_select           = "random",
                    output_path           = NULL,
                    hspaceTechniqueChoose = "gamma_dominant",
                    gamma_ratio_min       = 10,
                    constraint_threshold  = 0.1,
                    out_Decon             = NULL,
                    normalize             = TRUE,
                    seed                  = NULL) {

  # ---- Validation ----------------------------------------------------------
  args                  <- .validate_inputs(
    normalize             = normalize,
    algo_select           = algo_select,
    hspaceTechniqueChoose = hspaceTechniqueChoose,
    gamma_ratio_min       = gamma_ratio_min,
    cibersort_perm        = cibersort_perm,
    cibersort_QN          = cibersort_QN
  )
  normalize             <- args$normalize
  algo_select           <- args$algo_select
  hspaceTechniqueChoose <- args$hspaceTechniqueChoose
  gamma_ratio_min       <- args$gamma_ratio_min
  cibersort_perm        <- args$cibersort_perm
  cibersort_QN          <- args$cibersort_QN

  # ---- Preprocessing -------------------------------------------------------
  data      <- .prepare_data(reference, bulk, normalize)
  reference <- data$reference
  bulk      <- data$bulk

  # ---- Supervised deconvolution --------------------------------------------
  out_Dec <- .get_deconvolution(
    bulk             = bulk,
    reference        = reference,
    methodDeconv     = methodDeconv,
    cibersortx_email = cibersortx_email,
    cibersortx_token = cibersortx_token,
    out_Decon        = out_Decon,
    cibersort_perm   = cibersort_perm,
    cibersort_QN     = cibersort_QN
  )

  # ---- Data-set for NMF -----------------------------------------------------
  dataset <- list(B = bulk, W = reference, P = out_Dec)

  # ---- Output directory ----------------------------------------------------
  if (is.null(output_path)) output_path <- getwd()
  output_dir <- file.path(
    output_path,
    paste0("dicepro_", bulkName, "_", refName)
  )

  hp_params <- switch(
    hspaceTechniqueChoose,
    all                = c("gamma", "lambda_", "p_prime"),
    gamma_dominant = c("gamma", "lambda_", "p_prime"),
    restrictionEspace  = c("gamma", "lambda_factor", "p_prime")
  )

  seed <- seed %||% 42L

  # ---- hyper-parameter optimization -----------------------------------------
  out <- .run_hyperopt(
    dataset               = dataset,
    W_prime               = W_prime,
    bulkName              = bulkName,
    refName               = refName,
    hp_max_evals          = hp_max_evals,
    algo_select           = algo_select,
    output_base_dir       = output_path,
    hspaceTechniqueChoose = hspaceTechniqueChoose,
    output_dir            = output_dir,
    gamma_ratio_min       = gamma_ratio_min,
    constraint_threshold  = constraint_threshold,
    seed                  = seed
  )

  if (is.null(out)) {
    warning("No valid hyper-parameter configuration found.", call. = FALSE)
    return(invisible(NULL))
  }

  # ---- Save plots ----------------------------------------------------------
  out <- .save_outputs(out, output_dir, hp_params)

  structure(out, class = "dicepro")
}
