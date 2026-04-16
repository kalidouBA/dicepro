# =============================================================================
# dicepro Main Function
#
# Public API  : dicepro()
# Private fns :  .validate_inputs(), .normalize_zscore_per_gene()
# =============================================================================


# -----------------------------------------------------------------------------
# Internal input validation
# -----------------------------------------------------------------------------

#' Validate dicepro inputs
#'
#' Internal helper that validates and normalizes user parameters.
#'
#' @param normalize Logical. Whether to apply z-score normalization.
#' @param algo_select Character. Hyperparameter search strategy.
#' @param hspaceTechniqueChoose Character. Hyperparameter space strategy.
#'
#' @return A named list with validated parameters.
#'
#' @keywords internal
#' @noRd
.validate_inputs <- function(normalize,
                             algo_select,
                             hspaceTechniqueChoose) {

  if (!is.logical(normalize) || length(normalize) != 1L) {
    stop("'normalize' must be TRUE or FALSE.", call. = FALSE)
  }

  algo_select <- match.arg(
    tolower(algo_select),
    c("random", "tpe", "atpe", "anneal")
  )

  hspaceTechniqueChoose <- match.arg(
    hspaceTechniqueChoose,
    c("restrictionEspace", "all")
  )

  list(
    normalize = normalize,
    algo_select = algo_select,
    hspaceTechniqueChoose = hspaceTechniqueChoose
  )
}

# -----------------------------------------------------------------------------
# .normalize_zscore_per_gene  [private]
# -----------------------------------------------------------------------------

#' Z-score normalisation per gene
#'
#' Normalises each row (gene) of a numeric matrix by centering on its mean
#' and scaling by its standard deviation. Genes with SD = 0 or \code{NA}
#' are silently removed.
#'
#' @param mat A numeric matrix or data frame with genes as rows and samples
#'   as columns.
#'
#' @return A numeric matrix with the same columns as \code{mat} but
#'   potentially fewer rows (genes with SD = 0 removed). Each remaining
#'   row has mean = 0 and SD = 1.
#'
#' @importFrom stats sd
#' @keywords internal
#' @noRd
.normalize_zscore_per_gene <- function(mat) {

  mat <- as.matrix(mat)

  gene_sd    <- apply(mat, 1L, sd, na.rm = TRUE)
  keep_genes <- !is.na(gene_sd) & gene_sd > 0
  mat        <- mat[keep_genes, , drop = FALSE]

  t(apply(mat, 1L, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
}


# -----------------------------------------------------------------------------
# Internal input validation
# -----------------------------------------------------------------------------

#' Validate dicepro inputs
#'
#' Internal helper that validates and normalizes user parameters.
#'
#' @param normalize Logical. Whether to apply z-score normalization.
#' @param algo_select Character. Hyperparameter search strategy.
#' @param hspaceTechniqueChoose Character. Hyperparameter space strategy.
#'
#' @return A named list with validated parameters.
#'
#' @keywords internal
#' @noRd
.validate_inputs <- function(normalize,
                             algo_select,
                             hspaceTechniqueChoose) {

  if (!is.logical(normalize) || length(normalize) != 1L) {
    stop("'normalize' must be TRUE or FALSE.", call. = FALSE)
  }

  algo_select <- match.arg(
    tolower(algo_select),
    c("random", "tpe", "atpe", "anneal")
  )

  hspaceTechniqueChoose <- match.arg(
    hspaceTechniqueChoose,
    c("restrictionEspace", "all")
  )

  list(
    normalize = normalize,
    algo_select = algo_select,
    hspaceTechniqueChoose = hspaceTechniqueChoose
  )
}



# -----------------------------------------------------------------------------
# Data preprocessing
# -----------------------------------------------------------------------------

#' Prepare bulk and reference matrices
#'
#' Applies optional normalization and keeps intersecting genes only.
#'
#' @param reference Numeric matrix (genes x cell types).
#' @param bulk Numeric matrix (genes x samples).
#' @param normalize Logical.
#'
#' @return A list containing:
#' \describe{
#'   \item{reference}{Filtered reference matrix}
#'   \item{bulk}{Filtered bulk matrix}
#' }
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

  if (length(geneIntersect) == 0L) {
    stop("No common genes between 'reference' and 'bulk'.", call. = FALSE)
  }

  reference <- reference[geneIntersect, , drop = FALSE]
  bulk      <- bulk[geneIntersect, , drop = FALSE]

  rownames(reference) <- rownames(bulk) <- geneIntersect

  message(sprintf("Gene intersection: %d genes retained.", length(geneIntersect)))

  list(reference = reference, bulk = bulk)
}


# -----------------------------------------------------------------------------
# Deconvolution engine (private)
# -----------------------------------------------------------------------------

#' Run or validate supervised deconvolution
#'
#' Either executes \code{running_method()} or validates a user-provided matrix.
#'
#' @param bulk,reference Numeric matrices.
#' @param methodDeconv Character. Deconvolution method.
#' @param cibersortx_email,cibersortx_token Credentials.
#' @param out_Decon Optional precomputed matrix.
#'
#' @return Numeric matrix (cell types x samples).
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
      error = function(e) {
        stop(sprintf(
          "Invalid 'methodDeconv': '%s'. Must be one of: %s.",
          methodDeconv,
          paste(valid_methods, collapse = ", ")
        ), call. = FALSE)
      }
    )

    if (methodDeconv == "CSx" &&
        (is.null(cibersortx_email) || is.null(cibersortx_token))) {
      stop("CIBERSORTx credentials required for CSx.", call. = FALSE)
    }

    out_Decon <- running_method(
      bulk, reference,
      methodDeconv,
      cibersortx_email,
      cibersortx_token
    )

  } else {

    out_Decon <- as.matrix(out_Decon)

    if (ncol(out_Decon) != ncol(reference)) {
      stop("Mismatch: cell types vs reference.", call. = FALSE)
    }

    if (nrow(out_Decon) != ncol(bulk)) {
      stop("Mismatch: samples vs bulk.", call. = FALSE)
    }
  }

  t(out_Decon)
}


# -----------------------------------------------------------------------------
# Hyperparameter optimization (internal)
# -----------------------------------------------------------------------------

#' Run hyperparameter optimization pipeline
#'
#' Internal wrapper that executes the full hyperparameter search via
#' \code{run_experiment()} and selects the best configuration using
#' Pareto-front analysis through \code{best_hyperParams()}.
#'
#' @details
#' This function is part of the internal dicepro pipeline. It:
#' \enumerate{
#'   \item Runs hyperparameter exploration (\code{run_experiment})
#'   \item Collects all evaluated configurations (trials)
#'   \item Computes Pareto-optimal solutions
#'   \item Selects the best configuration via knee-point detection
#' }
#'
#' The selection criterion is handled internally by
#' \code{best_hyperParams()}.
#'
#' @param dataset Named list containing:
#'   \describe{
#'     \item{B}{Bulk expression matrix}
#'     \item{W}{Reference expression matrix}
#'     \item{P}{Initial deconvolution matrix}
#'   }
#' @param W_prime Numeric matrix or scalar.
#' Initial signature matrix for unknown cell types.
#' @param bulkName Character. Label for bulk dataset.
#' @param refName Character. Label for reference dataset.
#' @param hp_max_evals Integer. Number of hyperparameter configurations to evaluate.
#' @param algo_select Character. Sampling strategy for hyperparameter search:
#'   \code{"random"}, \code{"tpe"}, \code{"atpe"}, \code{"anneal"}.
#' @param output_base_dir Character. Root directory for all outputs.
#' @param hspaceTechniqueChoose Character. Hyperparameter space definition:
#'   \code{"all"} or \code{"restrictionEspace"}.
#' @param output_dir Character. Directory where final reports are saved.
#'
#' @return A list containing the best model configuration:
#' \describe{
#'   \item{hyperparameters}{Selected hyperparameters (lambda, gamma, etc.)}
#'   \item{metrics}{Loss and constraint values of the best solution}
#'   \item{trials}{Full hyperparameter search results}
#'   \item{W}{Estimated signature matrix}
#'   \item{H}{Estimated proportion matrix}
#'   \item{plot}{Pareto frontier visualization}
#' }
#'
#' @seealso
#' \code{\link{run_experiment}}, \code{\link{best_hyperParams}}
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
# Output persistence (internal)
# -----------------------------------------------------------------------------

#' Save hyperparameter optimization outputs
#'
#' Internal helper that generates and saves diagnostic plots from the
#' hyperparameter optimization step.
#'
#' @details
#' This function performs the following operations:
#' \enumerate{
#'   \item Creates a reporting directory under \code{output_dir/report}
#'   \item Generates a hyperparameter exploration plot using
#'     \code{plot_hyperopt()}
#'   \item Saves the hyperparameter optimization plot to PDF
#'   \item Saves the Pareto frontier plot to PDF
#'   \item Attaches the generated plot to the output object
#' }
#'
#' The function has side effects (file writing to disk).
#'
#' @param out List. Output object returned by the optimization pipeline.
#' Must contain at least:
#' \describe{
#'   \item{plot}{Pareto frontier ggplot object}
#' }
#'
#' @param output_dir Character. Root directory where results are stored.
#' A subdirectory \code{report/} will be created if it does not exist.
#'
#' @param hp_params Character vector. Names of hyperparameters used for
#' visualization (e.g. \code{c("gamma", "lambda_", "p_prime")}).
#'
#' @return The modified \code{out} object with an additional element:
#' \describe{
#'   \item{plot_hyperopt}{ggplot object of hyperparameter exploration}
#' }
#'
#' @seealso
#' \code{\link{plot_hyperopt}}
#'
#' @keywords internal
#' @noRd
.save_outputs <- function(out, output_dir, hp_params) {

  report_dir <- file.path(output_dir, "report")

  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
  }

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

#' Semi-supervised bulk RNA-seq deconvolution with hyperparameter optimization
#'
#' @description
#' Performs semi-supervised deconvolution of bulk RNA-seq data by combining:
#' a supervised estimation of known cell types and an unsupervised discovery
#' of latent components, with automatic hyperparameter optimization.
#'
#' @details
#' If \code{out_Decon} is provided, the supervised deconvolution step is skipped.
#' Otherwise, the method specified in \code{methodDeconv} is executed.
#'
#' Gene expression matrices are optionally normalized using z-score scaling
#' per gene, and only intersecting genes are retained.
#'
#' Hyperparameters are optimized using a Pareto-frontier-based strategy.
#'
#' @param reference Numeric matrix (genes × cell types)
#' @param bulk Numeric matrix (genes × samples)
#' @param methodDeconv Character. One of:
#'   \code{"CSx"}, \code{"DCQ"}, \code{"CDSeq"}, \code{"FARDEEP"}
#' @param cibersortx_email CIBERSORTx email (required for CSx)
#' @param cibersortx_token CIBERSORTx token (required for CSx)
#' @param W_prime Initial unknown signature matrix or 0
#' @param bulkName Label for bulk dataset
#' @param refName Label for reference dataset
#' @param hp_max_evals Number of hyperparameter evaluations
#' @param N_unknownCT Number of unknown cell types
#' @param algo_select Hyperparameter search strategy:
#'   \code{"random"}, \code{"tpe"}, \code{"atpe"}, \code{"anneal"}
#' @param output_path Output directory
#' @param hspaceTechniqueChoose Hyperparameter space:
#'   \code{"all"} or \code{"restrictionEspace"}
#' @param out_Decon Optional precomputed deconvolution matrix
#' @param normalize Logical. Apply z-score normalization per gene
#'
#' @return A list of class \code{"dicepro"} containing:
#' \itemize{
#'   \item \code{hyperparameters} selected lambda and gamma
#'   \item \code{metrics} loss and constraint values
#'   \item \code{trials} all evaluated configurations
#'   \item \code{W} estimated signatures
#'   \item \code{H} estimated proportions
#'   \item \code{plot} Pareto frontier plot
#'   \item \code{plot_hyperopt} hyperparameter space visualization
#' }
#'
#' @seealso
#' \code{\link{running_method}}, \code{\link{run_experiment}},
#' \code{\link{best_hyperParams}}
#'
#' @examples
#' \dontrun{
#' res <- dicepro(reference = ref_mat,
#'                 bulk = bulk_mat,
#'                 methodDeconv = "CSx",
#'                 cibersortx_email = "xxx",
#'                 cibersortx_token = "xxx")
#' }
#'
#' @export
dicepro <- function(reference, bulk,
                    methodDeconv = "CSx",
                    cibersortx_email = NULL,
                    cibersortx_token = NULL,
                    W_prime = 0,
                    bulkName = "Bulk",
                    refName = "Reference",
                    hp_max_evals = 100,
                    N_unknownCT = 1,
                    algo_select = "random",
                    output_path = NULL,
                    hspaceTechniqueChoose = "all",
                    out_Decon = NULL,
                    normalize = TRUE) {

  # ---- validation ----
  args <- .validate_inputs(normalize, algo_select, hspaceTechniqueChoose)
  normalize             <- args$normalize
  algo_select           <- args$algo_select
  hspaceTechniqueChoose <- args$hspaceTechniqueChoose

  # ---- preprocessing ----
  data <- .prepare_data(reference, bulk, normalize)
  reference <- data$reference
  bulk      <- data$bulk

  # ---- deconvolution ----
  out_Dec <- .get_deconvolution(
    bulk, reference,
    methodDeconv,
    cibersortx_email,
    cibersortx_token,
    out_Decon
  )

  # ---- dataset ----
  dataset <- list(B = bulk, W = reference, P = out_Dec)

  # ---- output dir ----
  dirName <- paste0("dicepro_", bulkName, "_", refName)
  if (is.null(output_path)) output_path <- getwd()
  output_dir <- file.path(output_path, dirName)

  hp_params <- switch(
    hspaceTechniqueChoose,
    all               = c("gamma", "lambda_", "p_prime"),
    restrictionEspace = c("gamma", "lambda_factor", "p_prime")
  )

  # ---- hyperopt ----
  out <- .run_hyperopt(
    dataset, W_prime,
    bulkName, refName,
    hp_max_evals,
    algo_select,
    output_path,
    hspaceTechniqueChoose,
    output_dir
  )

  if (is.null(out)) {
    warning("No valid hyperparameter configuration found.", call. = FALSE)
    return(invisible(NULL))
  }

  # ---- save ----
  out <- .save_outputs(out, output_dir, hp_params)

  structure(out, class = "dicepro")
}
