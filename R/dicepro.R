# =============================================================================
# dicepro.R — Main entry point
#
# Public API  : dicepro()
# Private fns : .validate_inputs(), .normalize_zscore_per_gene(),
#               .prepare_data(), .get_deconvolution(),
#               .run_hyperopt(), .save_outputs()
# =============================================================================


# -----------------------------------------------------------------------------
# .validate_inputs  [private]
# -----------------------------------------------------------------------------

#' Validate and normalise dicepro user parameters
#'
#' @param normalize           Logical.
#' @param algo_select         Character. Hyperparameter search strategy.
#' @param hspaceTechniqueChoose Character. Search-space strategy.
#'
#' @return Named list with validated \code{normalize}, \code{algo_select},
#'   and \code{hspaceTechniqueChoose}.
#'
#' @keywords internal
#' @noRd
.validate_inputs <- function(normalize,
                             algo_select,
                             hspaceTechniqueChoose) {

  if (!is.logical(normalize) || length(normalize) != 1L)
    stop("'normalize' must be TRUE or FALSE.", call. = FALSE)

  algo_select <- match.arg(
    tolower(algo_select),
    c("random", "tpe", "atpe", "anneal")
  )

  hspaceTechniqueChoose <- match.arg(
    hspaceTechniqueChoose,
    c("all", "restrictionEspace")
  )

  list(
    normalize             = normalize,
    algo_select           = algo_select,
    hspaceTechniqueChoose = hspaceTechniqueChoose
  )
}


# -----------------------------------------------------------------------------
# .normalize_zscore_per_gene  [private]
# -----------------------------------------------------------------------------

#' Z-score normalisation per gene
#'
#' Centres and scales each row (gene) of a numeric matrix. Genes with
#' \code{SD = 0} or \code{NA} are silently removed.
#'
#' @param mat Numeric matrix (genes × samples).
#' @return Numeric matrix with mean = 0 and SD = 1 per row.
#'
#' @keywords internal
#' @noRd
.normalize_zscore_per_gene <- function(mat) {

  mat       <- as.matrix(mat)
  gene_sd   <- apply(mat, 1L, stats::sd, na.rm = TRUE)
  keep      <- !is.na(gene_sd) & gene_sd > 0
  mat       <- mat[keep, , drop = FALSE]

  t(apply(mat, 1L, function(x) {
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  }))
}


# -----------------------------------------------------------------------------
# .prepare_data  [private]
# -----------------------------------------------------------------------------

#' Normalise and intersect bulk and reference matrices
#'
#' @param reference Numeric matrix (genes × cell types).
#' @param bulk      Numeric matrix (genes × samples).
#' @param normalize Logical.
#'
#' @return Named list: \code{reference}, \code{bulk} — both filtered to
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
#' @param bulk,reference       Numeric matrices (genes × ...).
#' @param methodDeconv         Character. Deconvolution method.
#' @param cibersortx_email,cibersortx_token CIBERSORTx credentials.
#' @param out_Decon            Optional precomputed matrix
#'   (samples × cell types) or \code{NULL}.
#'
#' @return Numeric matrix (cell types × samples).
#'
#' @keywords internal
#' @noRd
.get_deconvolution <- function(bulk,
                               reference,
                               methodDeconv,
                               cibersortx_email,
                               cibersortx_token,
                               out_Decon) {

  valid_methods <- c("CSx", "DCQ", "CDSeq", "FARDEEP")

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
      stop("CIBERSORTx credentials required for CSx.", call. = FALSE)

    out_Decon <- running_method(
      bulk, reference,
      methodDeconv,
      cibersortx_email,
      cibersortx_token
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

#' Execute the full hyperparameter optimisation pipeline
#'
#' Calls \code{\link{run_experiment}} then \code{\link{best_hyperParams}}.
#'
#' @param dataset             Named list: \code{B}, \code{W}, \code{P}.
#' @param W_prime             Numeric matrix or scalar.
#' @param bulkName,refName    Character scalars.
#' @param hp_max_evals        Positive integer.
#' @param algo_select         Character scalar.
#' @param output_base_dir     Character scalar.
#' @param hspaceTechniqueChoose Character scalar.
#' @param output_dir          Character scalar.
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
                          output_dir) {

  res <- run_experiment(
    dataset               = dataset,
    W_prime               = W_prime,
    bulkName              = bulkName,
    refName               = refName,
    hp_max_evals          = hp_max_evals,
    algo_select           = algo_select,
    output_base_dir       = output_base_dir,
    hspaceTechniqueChoose = hspaceTechniqueChoose
  )

  best_hyperParams(
    trials_df = res$trials,
    W         = res$W,
    H         = res$H,
    savePaths = output_dir
  )
}


# -----------------------------------------------------------------------------
# .save_outputs  [private]
# -----------------------------------------------------------------------------

#' Save hyperparameter optimisation diagnostic plots
#'
#' Creates \code{output_dir/report/}, saves two PDFs (hyperopt report and
#' Pareto frontier), and attaches the hyperopt plot to \code{out}.
#'
#' @param out        Named list returned by \code{.run_hyperopt}.
#' @param output_dir Character scalar. Root results directory.
#' @param hp_params  Character vector of hyperparameter names to plot.
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

#' Semi-supervised bulk RNA-seq deconvolution with hyperparameter optimisation
#'
#' @description
#' Combines supervised estimation of known cell types with unsupervised
#' discovery of latent components, with automatic Pareto-frontier-based
#' hyperparameter optimisation.
#'
#' @details
#' When \code{out_Decon} is provided the supervised step is skipped.
#' Gene matrices are optionally z-score normalised per gene, and only
#' intersecting genes are retained before optimisation.
#'
#' @param reference           Numeric matrix (genes × cell types).
#' @param bulk                Numeric matrix (genes × samples).
#' @param methodDeconv        Character. One of \code{"CSx"}, \code{"DCQ"},
#'   \code{"CDSeq"}, \code{"FARDEEP"}.
#' @param cibersortx_email    CIBERSORTx email (required for \code{"CSx"}).
#' @param cibersortx_token    CIBERSORTx token (required for \code{"CSx"}).
#' @param W_prime             Initial unknown-signature matrix or \code{0}.
#' @param bulkName            Character scalar. Label for the bulk dataset.
#' @param refName             Character scalar. Label for the reference.
#' @param hp_max_evals        Positive integer. Number of hyperparameter trials.
#' @param N_unknownCT         Positive integer. Number of unknown cell types.
#' @param algo_select         Character. One of \code{"random"}, \code{"tpe"},
#'   \code{"atpe"}, \code{"anneal"}.
#' @param output_path         Character scalar. Root output directory
#'   (\code{getwd()} by default).
#' @param hspaceTechniqueChoose Character. \code{"all"} or
#'   \code{"restrictionEspace"}.
#' @param out_Decon           Optional precomputed deconvolution matrix
#'   (samples × cell types).
#' @param normalize           Logical. Apply z-score normalisation per gene.
#'
#' @return An object of class \code{"dicepro"} (a named list) containing:
#' \itemize{
#'   \item \code{hyperparameters} — selected \eqn{\lambda} and \eqn{\gamma}
#'   \item \code{metrics}         — loss and constraint of the best trial
#'   \item \code{trials}          — all evaluated configurations
#'   \item \code{W}               — estimated unknown signature matrix
#'   \item \code{H}               — estimated proportion matrix
#'   \item \code{plot}            — Pareto frontier ggplot2 figure
#'   \item \code{plot_hyperopt}   — hyperparameter space ggplot2 figure
#' }
#' Returns \code{invisible(NULL)} with a warning when no valid configuration
#' is found.
#'
#' @seealso \code{\link{running_method}}, \code{\link{run_experiment}},
#'   \code{\link{best_hyperParams}}
#'
#' @examples
#' \dontrun{
#' res <- dicepro(
#'   reference    = BlueCode,
#'   bulk         = CellMixtures,
#'   methodDeconv = "FARDEEP",
#'   hp_max_evals = 50L
#' )
#' }
#'
#' @export
dicepro <- function(reference,
                    bulk,
                    methodDeconv          = "CSx",
                    cibersortx_email      = NULL,
                    cibersortx_token      = NULL,
                    W_prime               = 0,
                    bulkName              = "Bulk",
                    refName               = "Reference",
                    hp_max_evals          = 100L,
                    N_unknownCT           = 1L,
                    algo_select           = "random",
                    output_path           = NULL,
                    hspaceTechniqueChoose = "all",
                    out_Decon             = NULL,
                    normalize             = TRUE) {

  # ---- Validation ----------------------------------------------------------
  args                  <- .validate_inputs(normalize, algo_select, hspaceTechniqueChoose)
  normalize             <- args$normalize
  algo_select           <- args$algo_select
  hspaceTechniqueChoose <- args$hspaceTechniqueChoose

  # ---- Preprocessing -------------------------------------------------------
  data      <- .prepare_data(reference, bulk, normalize)
  reference <- data$reference
  bulk      <- data$bulk

  # ---- Supervised deconvolution --------------------------------------------
  out_Dec <- .get_deconvolution(
    bulk, reference,
    methodDeconv,
    cibersortx_email,
    cibersortx_token,
    out_Decon
  )

  # ---- Dataset for NMF -----------------------------------------------------
  dataset <- list(B = bulk, W = reference, P = out_Dec)

  # ---- Output directory ----------------------------------------------------
  if (is.null(output_path)) output_path <- getwd()
  output_dir <- file.path(
    output_path,
    paste0("dicepro_", bulkName, "_", refName)
  )

  hp_params <- switch(
    hspaceTechniqueChoose,
    all               = c("gamma", "lambda_", "p_prime"),
    restrictionEspace = c("gamma", "lambda_factor", "p_prime")
  )

  # ---- Hyperparameter optimisation -----------------------------------------
  out <- .run_hyperopt(
    dataset               = dataset,
    W_prime               = W_prime,
    bulkName              = bulkName,
    refName               = refName,
    hp_max_evals          = hp_max_evals,
    algo_select           = algo_select,
    output_base_dir       = output_path,
    hspaceTechniqueChoose = hspaceTechniqueChoose,
    output_dir            = output_dir
  )

  if (is.null(out)) {
    warning("No valid hyperparameter configuration found.", call. = FALSE)
    return(invisible(NULL))
  }

  # ---- Save plots ----------------------------------------------------------
  out <- .save_outputs(out, output_dir, hp_params)

  structure(out, class = "dicepro")
}
