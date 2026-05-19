#' Run CIBERSORTx Deconvolution Method
#'
#' This function runs the CIBERSORTx deconvolution method using Docker with the provided email and token.
#' The method is used to estimate cell type proportions in bulk RNA-seq data.
#'
#' @param bulk A matrix of bulk RNA-seq data.
#' @param reference A matrix of reference cell type gene expression profiles.
#' @param cibersortx_email The CIBERSORTx account email.
#' @param cibersortx_token The CIBERSORTx account token.
#' @param work_dir Character. Path to the working directory where CIBERSORTx
#'   output files are written. Defaults to \code{tempdir()}.
#'
#' @return A matrix containing the estimated cell type proportions.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Runs CIBERSORTx with the provided email and token.
#'   \item Reads the output results from CIBERSORTx.
#'   \item Extracts cell type proportions and returns them as a matrix.
#'}
#' @export

run_CSx <- function(bulk, reference,
                    cibersortx_email, cibersortx_token,
                    work_dir = tempdir()) {
  work_dir <- file.path(work_dir, paste0("dicepro_csx_", Sys.getpid()))
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  sig_file  <- file.path(work_dir, "SignCSx.txt")
  bulk_file <- file.path(work_dir, "bulk.txt")
  out_file  <- file.path(work_dir, "CIBERSORTx_Results.txt")

  sc_ref <- cbind.data.frame(GeneSymbol = rownames(reference), reference)
  utils::write.table(sc_ref, file = sig_file, sep = "\t",
                     row.names = FALSE, quote = FALSE)

  bulk_ <- cbind.data.frame(Gene = rownames(bulk), bulk)
  utils::write.table(bulk_, file = bulk_file, sep = "\t",
                     row.names = FALSE, quote = FALSE)

  command_to_run <- paste0(
    "docker run -v ", work_dir, ":/src/data ",
    "-v ", work_dir, ":/src/outdir cibersortx/fractions ",
    "--username ", cibersortx_email,
    " --token ",   cibersortx_token,
    " --mixture bulk.txt --sigmatrix SignCSx.txt --QN FALSE"
  )
  code <- system(command_to_run)

  if (code != 0L)
    stop("CIBERSORTx Docker run failed with exit code ", code, call. = FALSE)

  out_CSx   <- utils::read.delim2(out_file, row.names = 1)
  rownames_ <- rownames(out_CSx)
  out_CSx   <- as.data.frame(apply(
    out_CSx[, !names(out_CSx) %in% c("P.value", "Correlation", "RMSE")],
    2, as.numeric
  ))
  rownames(out_CSx) <- rownames_

  return(t(out_CSx))
}
