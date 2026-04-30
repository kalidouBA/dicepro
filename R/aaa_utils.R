# =============================================================================
# utils.R  — package-wide helpers, operators, and global variable declarations
#
# Private fns : convert_matrix_to_df(), .generate_experiment_paths(), `%||%`
# No public exports (all helpers are internal)
#
# NOTE: This file is named utils.R so R loads it early by convention.
#       All symbols defined here are available to every other file.
# =============================================================================


# -----------------------------------------------------------------------------
# Global variable declarations (R CMD CHECK NSE silencing)
# -----------------------------------------------------------------------------

# ggplot2 / dplyr column-name symbols used inside aes() or dplyr verbs
utils::globalVariables(c(
  "lambda_", "gamma", "frobNorm", "constraint", "abs_constraint",
  "p_prime", "duration", "penalty", "objectiveValue",
  "bound_label",          # lambda_gamma_initialize.R  → aes(colour = bound_label)
  "frob_H", "loss"        # hyperopt_selection.R → KraljicMatrix::get_frontier()
))


# -----------------------------------------------------------------------------
# convert_matrix_to_df  [private]
# -----------------------------------------------------------------------------

#' Convert a matrix to a data frame preserving dimnames
#'
#' @param mat A matrix. Non-matrix inputs are returned unchanged.
#' @return A data frame (same content, row names, column names) or \code{mat}
#'   unmodified.
#' @keywords internal
#' @noRd
convert_matrix_to_df <- function(mat) {
  if (!is.matrix(mat)) return(mat)
  df             <- as.data.frame(mat)
  rownames(df)   <- rownames(mat)
  colnames(df)   <- colnames(mat)
  df
}


# -----------------------------------------------------------------------------
# .generate_experiment_paths  [private]
# -----------------------------------------------------------------------------

#' Build and create standard experiment directory structure
#'
#' Creates \code{output_dir} and a sub-directory named
#' \code{dicepro_<bulkName>_<refName>} inside it, then returns both paths.
#'
#' @param output_dir Character scalar. Root output directory.
#' @param bulkName   Character scalar. Identifier for the bulk data-set.
#' @param refName    Character scalar. Identifier for the reference data-set.
#'
#' @return Named list with two elements:
#' \describe{
#'   \item{data_dir}{Path to the experiment subdirectory.}
#'   \item{config_path}{Path to the (not yet created) config JSON file.}
#' }
#'
#' @keywords internal
#' @noRd
.generate_experiment_paths <- function(output_dir, bulkName, refName) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  data_dir <- file.path(
    output_dir,
    paste0("dicepro_", bulkName, "_", refName)
  )
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  list(
    data_dir    = data_dir,
    config_path = file.path(data_dir, "optim.config.json")
  )
}


# -----------------------------------------------------------------------------
# %||%  [private]
# -----------------------------------------------------------------------------

#' Null-coalescing operator
#'
#' Returns \code{a} when it is not \code{NULL}, otherwise returns \code{b}.
#'
#' @param a Object to test.
#' @param b Default value returned when \code{a} is \code{NULL}.
#' @return \code{a} if \code{!is.null(a)}, else \code{b}.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
