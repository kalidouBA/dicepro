#' Run Cell Deconvolution
#'
#' Performs cell deconvolution on bulk RNA-seq data using multiple methods
#' such as CSx, DCQ, FARDEEP, and BayesPrism.
#'
#' @param bulk A numeric matrix of bulk RNA-seq data (rows = genes, columns = samples).
#' @param reference A numeric matrix of reference cell type gene expression profiles
#' (rows = genes, columns = cell types).
#' @param methodDeconv Character. The deconvolution method to use:
#' "CSx", "DCQ", "FARDEEP", or "BayesPrism".
#' @param cibersortx_email Optional. Email for CIBERSORTx account
#' (required if methodDeconv = "CSx").
#' @param cibersortx_token Optional. Token for CIBERSORTx account
#' (required if methodDeconv = "CSx").
#'
#' @return A numeric matrix with estimated cell type proportions for each sample.
#'
#' @details This function allows selection among different deconvolution algorithms.
#' The output is a matrix with normalized estimated cell fractions per sample.
#'
#' @seealso \code{run_CSx}, \code{dcq}, \code{fardeep}
#' @export
running_method <- function(
    bulk,
    reference,
    methodDeconv = "CSx",
    cibersortx_email = NULL,
    cibersortx_token = NULL
) {

  # Ensure numeric matrices while preserving dimnames
  dimname_ <- dimnames(reference)
  reference <- apply(reference, 2, as.numeric)
  dimnames(reference) <- dimname_

  dimname_ <- dimnames(bulk)
  bulk <- apply(bulk, 2, as.numeric)
  dimnames(bulk) <- dimname_

  # Keep only المشتر genes
  common <- intersect(row.names(bulk), row.names(reference))

  bulkIntersect <- as.data.frame(bulk[common, ])
  refIntersect  <- reference[common, ]

  row.names(bulkIntersect) <- common
  row.names(refIntersect)  <- common
  colnames(bulkIntersect)  <- colnames(bulk)

  # Run selected deconvolution method
  if (methodDeconv == "CSx") {

    out_Dec <- run_CSx(
      bulk,
      reference,
      cibersortx_email,
      cibersortx_token
    )

  } else if (methodDeconv == "DCQ") {

    out_Dec <- t(ComICS::dcq(
      reference_data   = refIntersect,
      mix_data         = bulkIntersect,
      marker_set       = as.data.frame(row.names(refIntersect)),
      alpha_used       = 0.05,
      lambda_min       = 0.2,
      number_of_repeats = 30
    )$average)

    out_Dec <- apply(out_Dec, 2, function(x) ifelse(x < 0, 0, x))

  } else if (methodDeconv == "FARDEEP") {

    out_Dec <- suppressWarnings(
      t(FARDEEP::fardeep(
        X         = reference,
        Y         = bulk,
        nn        = TRUE,
        intercept = TRUE,
        permn     = 100,
        QN        = FALSE
      )$abs.beta)
    )

  } else {

    stop("methodDeconv not supported")
  }

  # Normalize proportions
  out_Dec <- sweep(out_Dec, 2, colSums(out_Dec), FUN = "/")

  return(out_Dec)
}
