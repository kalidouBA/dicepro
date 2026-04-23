# =============================================================================
# simulation.R
# Core simulation functions for bulk RNA-seq deconvolution benchmarking.
#
# Public functions:
#   - generateProp()          : cell-type proportion matrix generation
#   - generate_ref_matrix()   : synthetic reference signature matrix
#   - simulation()            : bulk RNA-seq simulation with noise
#   - simulation_bluecode()   : bulk simulation using the BlueCode reference
#
# Internal helpers (not exported):
#   - .mvrnorm()              : multivariate normal sampler
#   - .rdirichlet()           : Dirichlet sampler
#   - .rcorrmatrix()          : random correlation matrix generator
#   - .bluecode_cell_groups   : canonical cell-group definition for BlueCode
#   - .generate_bluecode_prop(): hierarchical proportion sampler for BlueCode
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' @keywords internal
.mvrnorm <- function(n, mu, Sigma, empirical = FALSE) {
  p     <- length(mu)
  eS    <- eigen(Sigma, symmetric = TRUE)
  ev    <- eS$values
  if (!all(ev >= -1e-06 * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, scale = FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, scale = FALSE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  t(X)
}

#' @keywords internal
.rdirichlet <- function(n, alpha) {
  k     <- length(alpha)
  draws <- matrix(stats::rgamma(n * k, shape = alpha, rate = 1L),
                  nrow = n, ncol = k, byrow = TRUE)
  draws / rowSums(draws)
}

#' @keywords internal
.rcorrmatrix <- function(d) {
  S <- diag(d)
  P <- matrix(0, d, d)

  for (i in seq_len(d - 1)) {
    for (j in seq(i + 1, d)) {
      P[i, j] <- runif(1, -1, 1)
      p <- P[i, j]
      if (i > 1) {
        for (k in seq_len(i - 1)) {
          p <- p * sqrt((1 - P[k, j]^2) * (1 - P[k, i]^2)) + P[k, i] * P[k, j]
        }
      }
      S[i, j] <- p
      S[j, i] <- p
    }
  }

  eS <- eigen(S, symmetric = TRUE)
  ev <- pmax(eS$values, 1e-8)
  S  <- eS$vectors %*% diag(ev) %*% t(eS$vectors)
  D  <- diag(1 / sqrt(diag(S)))
  D %*% S %*% D
}


# -----------------------------------------------------------------------------
# BlueCode cell-group definition
#
# Maps the 34 cell types in colnames(BlueCode) to their five tissue
# compartments. This object is the single source of truth for the
# hierarchical structure used in simulation_bluecode() and in
# generateProp(scenario = "bluecode").
#
# Compartment membership and ordering must stay aligned with colnames(BlueCode).
# -----------------------------------------------------------------------------

#' @keywords internal
.bluecode_cell_groups <- list(
  Immune = c(
    "B.cell.naive", "B.cell.memory",
    "T.Cell.CD4", "T.cell.CD4.memory", "T.cell.CD8.memory",
    "NK.cell.CD56",
    "Monocyte.CD14", "Macrophage", "Mature.neutrophil"
  ),
  Stromal = c(
    "Fibroblast.cardiac.ventricle",
    "Fibroblast.of.arm",
    "Fibroblast.of.lung",
    "Fibroblast.of.the.aortic.adventitia",
    "Fibroblast.or.papilla.dermal.cell",
    "MSC.like.pluripotent.cell",
    "Chondrocyte.articular",
    "Osteoblast"
  ),
  Endothelial = c(
    "Endothelial.large.blood.vessel",
    "Endothelial.microvascular.mammary.or.endometrial",
    "Endothelial.microvascular.non.reproductive"
  ),
  Epithelial = c(
    "Epithelial.cell.mammary",
    "Epithelial.renal.cortical",
    "Epithelial.respiratory",
    "Hair.follicular.keratinocyte",
    "Melanocyte.of.skin"
  ),
  Muscle = c(
    "Smooth.muscle.cell.aortic",
    "Smooth.muscle.cell.bronchial",
    "Smooth.muscle.cell.coronary.artery",
    "Smooth.muscle.cell.pulmonary.artery",
    "Smooth.muscle.cell.of.bladder",
    "Smooth.muscle.cell.of.trachea",
    "Smooth.muscle.cell.uterine",
    "Myocyte.regular.cardiac",
    "Myometrial.cell"
  )
)


# -----------------------------------------------------------------------------
# Internal proportion sampler for the BlueCode hierarchical scenario
#
# Implements a two-level Dirichlet model:
#   Level 1 — tissue compartment proportions drawn from Dirichlet(alpha_groups)
#   Level 2 — cell-type sub-proportions within each compartment drawn from
#              Dirichlet(rep(alpha_subtypes[[grp]], n_subtypes)) and scaled
#              by the compartment proportion.
#
# Unlike the generic "hierarchical" scenario in generateProp(), this function:
#   - preserves real cell-type names from BlueCode (no CellType_k renaming)
#   - uses compartment-specific alpha scalars (alpha_subtypes list)
#   - guarantees column ordering consistent with .bluecode_cell_groups
# -----------------------------------------------------------------------------

#' @keywords internal
.generate_bluecode_prop <- function(nSample,
                                    cell_groups,
                                    alpha_groups,
                                    alpha_subtypes,
                                    seed) {
  set.seed(seed)

  # Level 1: draw compartment proportions for all samples simultaneously
  group_props <- .rdirichlet(nSample, alpha_groups)
  colnames(group_props) <- names(cell_groups)
  rownames(group_props) <- paste0("sample_", seq_len(nSample))

  # Pre-allocate the full proportion matrix (samples x cell types)
  all_celltypes <- unlist(cell_groups, use.names = FALSE)
  prop_matrix   <- matrix(
    0,
    nrow     = nSample,
    ncol     = length(all_celltypes),
    dimnames = list(
      paste0("sample_", seq_len(nSample)),
      all_celltypes
    )
  )

  # Level 2: for each compartment, draw sub-type proportions and scale
  for (grp in names(cell_groups)) {
    subtypes <- cell_groups[[grp]]
    alpha_g  <- alpha_subtypes[[grp]]

    # One independent Dirichlet draw per sample (nSample x n_subtypes)
    sub_props <- .rdirichlet(nSample, rep(alpha_g, length(subtypes)))
    colnames(sub_props) <- subtypes

    # Weight sub-type proportions by the compartment proportion
    prop_matrix[, subtypes] <- sub_props * group_props[, grp]
  }

  prop_matrix
}


# =============================================================================
# Public functions
# =============================================================================

#' Generate Cell-Type Proportion Matrix
#'
#' Simulates cell-type proportion matrices for bulk RNA-seq deconvolution
#' benchmarking. Four scenarios are supported: equal proportions across cell
#' types, uniform random sampling, hierarchical Dirichlet sampling reflecting
#' realistic tissue compartment organization, or a BlueCode-specific
#' hierarchical model that preserves real cell-type names.
#'
#' @param n_cell_types Integer. Number of distinct cell types. Ignored when
#'   \code{scenario = "bluecode"} (fixed at 34 by the BlueCode reference).
#' @param nSample      Integer. Number of samples to generate.
#' @param nCell        Integer. Total number of cells per sample (used for
#'   rounding precision in the uniform scenario).
#' @param scenario     Character. Proportion generation strategy. One of:
#'   \code{"even"} for near-equal proportions,
#'   \code{"uniform"} for random uniform sampling,
#'   \code{"hierarchical"} for two-level Dirichlet sampling with generic
#'   \code{CellType_k} names, or
#'   \code{"bluecode"} for two-level Dirichlet sampling using the real
#'   cell-type names and compartment structure of the BlueCode reference.
#'   Any other value defaults to a flat Dirichlet draw.
#' @param alpha_groups Named numeric vector of length 5. Dirichlet
#'   concentration parameters for the five tissue compartments
#'   (Immune, Stromal, Endothelial, Epithelial, Muscle). Only used when
#'   \code{scenario = "bluecode"}.
#'   Default: \code{c(Immune=4, Stromal=2.5, Endothelial=1.8, Epithelial=1.8,
#'   Muscle=1.5)}.
#' @param alpha_subtypes Named list with one scalar per compartment. Controls
#'   the concentration of sub-type Dirichlet draws within each compartment.
#'   Higher values produce more even sub-type distributions. Only used when
#'   \code{scenario = "bluecode"}.
#'   Default: all compartments set to \code{8}.
#' @param seed Integer. Random seed. Only used when \code{scenario =
#'   "bluecode"}. Default: \code{1234}.
#'
#' @return A numeric matrix of dimensions \code{nSample x n_cell_types} where
#'   each row sums to 1. For all scenarios except \code{"bluecode"}, row names
#'   are \code{Sample_1}, ..., \code{Sample_nSample} and column names are
#'   \code{CellType_1}, ..., \code{CellType_n_cell_types}. For
#'   \code{"bluecode"}, row names are \code{sample_1}, ...,
#'   \code{sample_nSample} and column names are the real BlueCode cell-type
#'   names ordered by compartment.
#'
#' @export
#' @importFrom stats rnorm runif aggregate
#'
#' @examples
#' if (interactive()) {
#'   set.seed(2101)
#'
#'   # Equal proportions
#'   prop_even <- generateProp(nSample = 20, n_cell_types = 10,
#'     scenario = "even")
#'
#'   # Generic hierarchical Dirichlet
#'   prop_hier <- generateProp(nSample = 20, n_cell_types = 34,
#'     scenario = "hierarchical")
#'   all(abs(rowSums(prop_hier) - 1) < 1e-8)  # TRUE
#'
#'   # BlueCode hierarchical Dirichlet (real cell-type names)
#'   prop_bc <- generateProp(nSample = 20, scenario = "bluecode", seed = 42)
#'   colnames(prop_bc)  # real BlueCode cell-type names
#' }
generateProp <- function(n_cell_types   = NULL,
                         nSample,
                         nCell          = 500,
                         scenario       = NULL,
                         alpha_groups   = c(
                           Immune      = 4.0,
                           Stromal     = 2.5,
                           Endothelial = 1.8,
                           Epithelial  = 1.8,
                           Muscle      = 1.5
                         ),
                         alpha_subtypes = list(
                           Immune      = 8,
                           Stromal     = 8,
                           Endothelial = 8,
                           Epithelial  = 8,
                           Muscle      = 8
                         ),
                         seed = 1234) {

  if (scenario == "even") {
    # ── Near-equal proportions with small Gaussian perturbation ──────────
    m <- round(
      matrix(
        abs(rnorm(n_cell_types,
                  mean = 1 / n_cell_types,
                  sd   = 0.01)),
        ncol = n_cell_types
      ), 3
    )

  } else if (scenario == "uniform") {
    # ── Uniform random sampling with variable number of active cell types ─
    min.percentage <- 1
    max.percentage <- 99
    CT             <- paste0("CellType_", seq(n_cell_types))
    P_all          <- data.frame()

    for (i in seq_len(nSample)) {
      num.CT      <- sample(x = 2:length(CT), 1)
      selected.CT <- sample(CT, num.CT, replace = FALSE)

      P <- runif(num.CT, min.percentage, max.percentage)
      P <- round(P / sum(P), digits = log10(nCell))
      P <- data.frame(CT       = selected.CT,
                      expected = P,
                      stringsAsFactors = FALSE)

      missing.CT <- data.frame(
        CT       = CT[!CT %in% selected.CT],
        expected = rep(0, length(CT) - num.CT),
        stringsAsFactors = FALSE
      )

      P_all <- rbind.data.frame(P_all, rbind(P, missing.CT))
    }

    m <- matrix(
      aggregate(P_all$expected,
                list(P_all$CT),
                FUN = mean)$x,
      ncol = n_cell_types
    )

  } else if (scenario == "hierarchical") {
    # ── Generic two-level hierarchical Dirichlet (CellType_k names) ───────
    # Level 1: major tissue compartments
    alpha_groups_hier <- c(
      Immune      = 6.0,
      Stromal     = 2.0,
      Endothelial = 1.5,
      Epithelial  = 1.5,
      Muscle      = 1.0
    )
    groups <- .rdirichlet(nSample, alpha_groups_hier)
    colnames(groups) <- names(alpha_groups_hier)

    # Level 2: compartment-specific sub-type alpha vectors
    alpha_immune <- c(
      B_cell_naive      = 1.8,
      B_cell_memory     = 1.2,
      T_CD4_naive       = 2.0,
      T_CD4_memory      = 1.5,
      T_CD8_memory      = 1.5,
      NK_cell           = 1.0,
      Monocyte_CD14     = 1.0,
      Macrophage        = 1.0,
      Mature_neutrophil = 1.0
    )
    immune      <- .rdirichlet(nSample, alpha_immune)        * groups[, "Immune"]
    stromal     <- .rdirichlet(nSample, rep(2, 8))           * groups[, "Stromal"]
    endothelial <- .rdirichlet(nSample, rep(2, 3))           * groups[, "Endothelial"]
    epithelial  <- .rdirichlet(nSample, c(2, 2, 2, 1.5, 1.5)) * groups[, "Epithelial"]
    muscle      <- .rdirichlet(nSample, rep(1, 9))           * groups[, "Muscle"]

    m <- cbind(immune, stromal, endothelial, epithelial, muscle)

    # Trim or pad to match the requested number of cell types
    if (ncol(m) > n_cell_types)
      m <- m[, seq_len(n_cell_types), drop = FALSE]

    if (ncol(m) < n_cell_types) {
      extra <- .rdirichlet(nSample, rep(1, n_cell_types - ncol(m)))
      m     <- cbind(m, extra)
    }
    n_cell_types <- ncol(m)

  } else if (scenario == "bluecode") {
    # ── BlueCode hierarchical Dirichlet (real cell-type names) ────────────
    # Delegates entirely to .generate_bluecode_prop(); returns early because
    # row/column naming is already handled there.
    m <- .generate_bluecode_prop(
      nSample        = nSample,
      cell_groups    = .bluecode_cell_groups,
      alpha_groups   = alpha_groups,
      alpha_subtypes = alpha_subtypes,
      seed           = seed
    )
    return(m)   # dimnames already set; skip the generic renaming below

  } else {
    # ── Flat Dirichlet draw (default fallback) ────────────────────────────
    m <- .rdirichlet(n = nSample, alpha = rep(1, n_cell_types))
  }

  # Generic row/column naming for all non-bluecode scenarios
  dimnames(m) <- list(
    paste0("Sample_",   seq_len(nrow(m))),
    paste0("CellType_", seq_len(ncol(m)))
  )

  return(m)
}


#' Generate a Reference Signature Matrix
#'
#' Constructs a synthetic gene-by-cell-type reference matrix using either
#' Poisson-based (via Gaussian copula) or log-normal gene expression models,
#' with optional block sparsity and TPM normalization.
#'
#' @param loi             Character. Expression model: \code{"rpois"} for
#'   Poisson-based simulation via Gaussian copula or any other value for
#'   log-normal simulation. Default: \code{"rpois"}.
#' @param tpm             Logical. If \code{TRUE}, columns are TPM-normalized.
#'   Default: \code{FALSE}.
#' @param bloc            Logical. If \code{TRUE}, introduces block sparsity
#'   so that each cell type expresses only a subset of genes.
#'   Default: \code{FALSE}.
#' @param nGenesByCellType Integer. Number of genes per cell-type block
#'   (used when \code{bloc = TRUE}).
#' @param nCell           Integer. Total number of cells (unused directly;
#'   retained for API compatibility).
#' @param nCellsType      Integer. Number of cell types. Default: \code{10}.
#' @param nGenes          Integer. Number of genes. Default: \code{500}.
#' @param lambda_vec      Numeric vector of length \code{nCellsType}. Per
#'   cell-type Poisson lambda parameters. If \code{NULL}, drawn uniformly in
#'   \code{[20, 60]}. Only used when \code{loi = "rpois"}.
#' @param corr            Numeric matrix. Inter-cell-type correlation matrix.
#'   If \code{NULL}, generated via \code{.rcorrmatrix}.
#' @param sparse          Logical. If \code{TRUE}, introduces random sparsity
#'   by zeroing entries with probability \code{1 - prob_sparse}.
#'   Default: \code{FALSE}.
#' @param prob_sparse     Numeric. Probability of a non-zero entry when
#'   \code{sparse = TRUE}. Default: \code{NULL}.
#'
#' @return A numeric matrix of dimensions \code{nGenes x nCellsType}.
#'   Row names are \code{Gene_1}, ..., \code{Gene_nGenes}; column names are
#'   \code{CellType_1}, ..., \code{CellType_nCellsType}.
#'
#' @details
#' For \code{loi = "rpois"}, correlated Poisson counts are generated via a
#' Gaussian copula: multivariate normal draws with covariance \code{corr} are
#' mapped to uniform marginals via \code{pnorm}, then to Poisson quantiles via
#' \code{qpois}. This preserves the inter-cell-type correlation structure
#' without requiring \pkg{SimMultiCorrData}.
#'
#' For any other value of \code{loi}, a log-normal model is used: multivariate
#' normal draws are exponentiated after a location shift of 6.
#'
#' @export
#' @importFrom stats runif rbinom rnorm pnorm qpois
#'
#' @examples
#' if (interactive()) {
#'   set.seed(2101)
#'   ref_pois  <- generate_ref_matrix(loi = "rpois", nGenes = 50, nCellsType = 5)
#'   ref_gauss <- generate_ref_matrix(loi = "gauss", nGenes = 50, nCellsType = 5)
#'   dim(ref_pois)   # 50 x 5
#'   dim(ref_gauss)  # 50 x 5
#' }
generate_ref_matrix <- function(loi              = "rpois",
                                tpm              = FALSE,
                                bloc             = FALSE,
                                nGenesByCellType = 50,
                                nCell            = 500,
                                nCellsType       = 10,
                                nGenes           = 500,
                                lambda_vec       = NULL,
                                corr             = NULL,
                                sparse           = FALSE,
                                prob_sparse      = NULL) {

  if (is.null(corr))
    corr <- .rcorrmatrix(nCellsType)

  if (loi == "rpois") {
    # Gaussian copula approach: MVN -> uniform marginals -> Poisson quantiles
    eS <- eigen(corr, symmetric = TRUE)
    ev <- pmax(eS$values, 0)
    Z  <- matrix(rnorm(nGenes * nCellsType), nrow = nGenes) %*%
      t(eS$vectors %*% diag(sqrt(ev), nCellsType))
    U <- pnorm(Z)

    if (is.null(lambda_vec))
      lambda_vec <- runif(nCellsType, min = 20, max = 60)

    if (length(lambda_vec) != nCellsType)
      stop("'lambda_vec' must have length equal to 'nCellsType'.")

    counts <- matrix(
      mapply(
        function(u, lam) qpois(u, lambda = lam),
        as.vector(U),
        rep(lambda_vec, each = nGenes)
      ),
      nrow = nGenes,
      ncol = nCellsType
    )

  } else {
    # Log-normal model: exponentiate shifted MVN draws
    mu     <- runif(nCellsType, -1, 0.5)
    counts <- exp(
      .mvrnorm(n = nGenes, mu = mu, Sigma = corr, empirical = TRUE) + 6
    )
  }

  colnames(counts) <- paste0("CellType_", seq_len(nCellsType))
  rownames(counts) <- paste0("Gene_",     seq_len(nGenes))

  if (sparse)
    counts <- apply(counts, 2, function(x)
      x * rbinom(n = nGenes, size = 1, prob = prob_sparse))

  if (bloc) {
    listNGE <- seq(1, nGenes + 1, nGenesByCellType)
    for (ind in seq_len(nCellsType))
      counts[-(listNGE[ind]:(listNGE[ind + 1] - 1)), ind] <- 0
  }

  if (tpm)
    counts <- t(1e6 * t(counts) / colSums(counts))

  return(as.matrix(counts))
}


#' Simulate Bulk RNA-seq Data with Biological and Technical Noise
#'
#' Generates synthetic bulk RNA-seq expression matrices by linearly mixing
#' reference cell-type profiles according to simulated proportions, with
#' optional biological and technical noise. The noise model follows the
#' Gaussian likelihood framework underlying the dicepro objective function:
#' biological noise is multiplicative, whereas technical noise is additive.
#' Consequently, normalization to a Gaussian scale prior
#' to deconvolution is required.
#'
#' @param W               Numeric matrix (genes x cell types) or \code{NULL}.
#'   If \code{NULL}, a reference matrix is generated internally via
#'   \code{\link{generate_ref_matrix}}.
#' @param prop            Numeric matrix (samples x cell types) or \code{NULL}.
#'   If \code{NULL}, proportions are generated via \code{\link{generateProp}}.
#' @param nSample         Integer. Number of samples. Default: \code{50}.
#' @param nCell           Integer. Total cells per sample. Default: \code{500}.
#' @param nCellsType      Integer. Number of cell types. Default: \code{50}.
#' @param nGenes          Integer. Number of genes. Default: \code{500}.
#' @param lambda_vec      Numeric vector of length \code{nCellsType}. Per
#'   cell-type Poisson lambda parameters. Passed to
#'   \code{\link{generate_ref_matrix}}.
#' @param corr            Numeric matrix. Inter-cell-type correlation.
#'   If \code{NULL}, generated internally.
#' @param scenario        Character. Proportion scenario passed to
#'   \code{\link{generateProp}}. Use \code{"hierarchical"} for realistic
#'   tissue composition (recommended).
#' @param loi             Character. Expression law for reference generation.
#'   \code{"rpois"} for Poisson via Gaussian copula, any other value for
#'   log-normal. Default: \code{"rpois"}.
#' @param tpm             Logical. TPM normalization of reference.
#'   Default: \code{FALSE}.
#' @param bloc            Logical. Block sparsity in reference.
#'   Default: \code{FALSE}.
#' @param nGenesByCellType Integer. Genes per block (when \code{bloc = TRUE}).
#'   Default: \code{50}.
#' @param sparse          Logical. Random sparsity in reference.
#'   Default: \code{FALSE}.
#' @param prob_sparse     Numeric. Sparsity probability. Default: \code{0.5}.
#' @param sigma_bio       Numeric. Standard deviation of multiplicative
#'   biological noise. Default: \code{0.07}.
#' @param sigma_tech      Numeric. Standard deviation of additive technical
#'   noise (proportional to expression level). Default: \code{0.07}.
#' @param seed            Integer. Random seed. Default: \code{1234}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{p}}{data.frame. True cell-type proportions
#'     (samples x cell types).}
#'   \item{\code{W}}{data.frame. Reference signature matrix
#'     (genes x cell types).}
#'   \item{\code{B}}{data.frame. Noisy bulk expression matrix
#'     (genes x samples).}
#' }
#'
#' @details
#' The bulk expression matrix is generated as:
#' \deqn{\mathbf{B} = \mathbf{W} \cdot \mathbf{P}^\top}
#' Noise is then applied as:
#' \deqn{B^*_{gn} = B_{gn}(1 + \varepsilon^{\text{bio}}_{gn}) +
#'       \varepsilon^{\text{tech}}_{gn} \cdot B_{gn}}
#' where \eqn{\varepsilon^{\text{bio}} \sim \mathcal{N}(0, \sigma^2_{\text{bio}})}
#' and \eqn{\varepsilon^{\text{tech}} \sim \mathcal{N}(0, \sigma^2_{\text{tech}})}.
#' Negative values are clipped to zero.
#'
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' if (interactive()) {
#'   set.seed(2101)
#'   sim <- simulation(
#'     scenario   = "hierarchical",
#'     nSample    = 20,
#'     nGenes     = 100,
#'     nCellsType = 10,
#'     sigma_bio  = 0.07,
#'     sigma_tech = 0.07,
#'     seed       = 2101
#'   )
#'   dim(sim$p)  # 20 x 10
#'   dim(sim$B)  # 100 x 20
#' }
simulation <- function(W                = NULL,
                       prop             = NULL,
                       nSample          = 50,
                       nCell            = 500,
                       nCellsType       = 50,
                       nGenes           = 500,
                       lambda_vec       = NULL,
                       corr             = NULL,
                       scenario         = NULL,
                       loi              = "rpois",
                       tpm              = FALSE,
                       bloc             = FALSE,
                       nGenesByCellType = 50,
                       sparse           = FALSE,
                       prob_sparse      = 0.5,
                       sigma_bio        = 0.07,
                       sigma_tech       = 0.07,
                       seed             = 1234) {

  # ── Reference matrix ──────────────────────────────────────────────────────
  if (is.null(W)) {
    W <- generate_ref_matrix(
      loi              = loi,
      nCellsType       = nCellsType,
      nGenes           = nGenes,
      lambda_vec       = lambda_vec,
      corr             = corr,
      tpm              = tpm,
      bloc             = bloc,
      nGenesByCellType = nGenesByCellType,
      sparse           = sparse,
      prob_sparse      = prob_sparse
    )
  } else {
    nGenes     <- nrow(W)
    nCellsType <- ncol(W)
    dimnames_  <- dimnames(W)
    W          <- apply(W, 2, as.numeric)
    dimnames(W) <- dimnames_
  }

  # ── Proportion matrix ─────────────────────────────────────────────────────
  if (is.null(prop)) {
    prop <- generateProp(
      n_cell_types = nCellsType,
      nSample      = nSample,
      nCell        = nCell,
      scenario     = scenario
    )
  }
  prop <- as.matrix(prop)

  # ── Linear mixing: B = W %*% t(prop) ─────────────────────────────────────
  bulk <- as.matrix(W) %*% t(prop)
  rownames(bulk) <- rownames(W)

  set.seed(seed)

  # ── Biological noise: multiplicative, mean-zero Gaussian ──────────────────
  bio_noise  <- matrix(
    rnorm(length(bulk), mean = 0, sd = sigma_bio),
    nrow = nrow(bulk), ncol = ncol(bulk)
  )
  bulk_noisy <- bulk * (1 + bio_noise)

  # ── Technical noise: additive, scaled to expression level ─────────────────
  tech_noise <- matrix(
    rnorm(length(bulk), mean = 0, sd = sigma_tech),
    nrow = nrow(bulk), ncol = ncol(bulk)
  )
  bulk_noisy <- bulk_noisy + tech_noise * bulk

  # ── Clip negative values ──────────────────────────────────────────────────
  bulk_noisy[bulk_noisy < 0] <- 0

  rownames(bulk_noisy) <- rownames(bulk)
  colnames(bulk_noisy) <- colnames(bulk)
  colnames(prop)       <- colnames(W)
  rownames(prop)       <- colnames(bulk_noisy) <- paste0("sample_", seq_len(nSample))

  return(list(
    p = prop,
    W = as.data.frame(W),
    B = as.data.frame(bulk_noisy)
  ))
}


#' Simulate Bulk RNA-seq Data Using the BlueCode Reference Matrix
#'
#' Generates synthetic bulk RNA-seq expression matrices by mixing the bundled
#' BlueCode reference signature matrix with hierarchically simulated cell-type
#' proportions, then applying biological and technical noise.
#'
#' The proportion model follows a two-level Dirichlet hierarchy that mirrors
#' the biological organization of human tissues:
#' \enumerate{
#'   \item \strong{Compartment level} — five major tissue compartments
#'     (Immune, Stromal, Endothelial, Epithelial, Muscle) are drawn from
#'     \eqn{\text{Dirichlet}(\alpha_{\text{groups}})}.
#'   \item \strong{Cell-type level} — within each compartment, sub-type
#'     proportions are drawn from a symmetric Dirichlet controlled by
#'     \code{alpha_subtypes} and scaled by the compartment proportion.
#' }
#' This two-level structure induces realistic positive within-compartment and
#' negative between-compartment correlations between cell types.
#'
#' @param nSample       Integer. Number of samples to simulate.
#'   Default: \code{50}.
#' @param alpha_groups  Named numeric vector of length 5. Dirichlet
#'   concentration parameters for the five tissue compartments
#'   (Immune, Stromal, Endothelial, Epithelial, Muscle). Larger values
#'   produce more even compartment proportions.
#'   Default: \code{c(Immune=4, Stromal=2.5, Endothelial=1.8,
#'   Epithelial=1.8, Muscle=1.5)}.
#' @param alpha_subtypes Named list with one numeric scalar per compartment.
#'   Controls the concentration of sub-type Dirichlet draws within each
#'   compartment. Higher values produce more even sub-type distributions.
#'   Default: all compartments set to \code{8}.
#' @param sigma_bio     Numeric. Standard deviation of multiplicative
#'   biological noise (log-normal). Default: \code{0.15}.
#' @param sigma_tech    Numeric. Standard deviation of additive technical
#'   noise, expressed as a fraction of the global expression standard
#'   deviation. Default: \code{0.02}.
#' @param seed          Integer. Random seed for full reproducibility.
#'   Proportion sampling uses \code{seed}; noise sampling uses
#'   \code{seed + 1L} so that both components can be varied independently.
#'   Default: \code{1234}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{p}}{data.frame (samples x 34 cell types). True cell-type
#'     proportions with real BlueCode cell-type names. Each row sums to 1.}
#'   \item{\code{W}}{data.frame (genes x 34 cell types). BlueCode reference
#'     signature matrix, column-ordered to match \code{p}.}
#'   \item{\code{B}}{data.frame (genes x samples). Noisy simulated bulk
#'     expression matrix.}
#' }
#'
#' @details
#' The bulk matrix is computed as \eqn{B = W \cdot P^\top}, then noise is
#' applied as:
#' \deqn{B^*_{gn} = B_{gn} \cdot \exp(\varepsilon^{\text{bio}}_{gn}) +
#'       \varepsilon^{\text{tech}}_{gn}}
#' where
#' \eqn{\varepsilon^{\text{bio}} \sim \mathcal{N}(0, \sigma^2_{\text{bio}})}
#' and
#' \eqn{\varepsilon^{\text{tech}} \sim
#'   \mathcal{N}(0,\, \sigma^2_{\text{tech}} \cdot \mathrm{sd}(B))}.
#' Negative values are clipped to zero.
#'
#' The function validates at runtime that the cell-type names in
#' \code{.bluecode_cell_groups} exactly match \code{colnames(BlueCode)},
#' ensuring consistency between the proportion model and the reference matrix.
#'
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' if (interactive()) {
#'   sim <- simulation_bluecode(
#'     nSample    = 30,
#'     sigma_bio  = 0.15,
#'     sigma_tech = 0.02,
#'     seed       = 42
#'   )
#'   dim(sim$p)                            # 30 x 34
#'   dim(sim$W)                            # nGenes x 34
#'   dim(sim$B)                            # nGenes x 30
#'   all(abs(rowSums(sim$p) - 1) < 1e-8)  # TRUE
#' }
simulation_bluecode <- function(
    nSample        = 50,
    alpha_groups   = c(
      Immune      = 4.0,
      Stromal     = 2.5,
      Endothelial = 1.8,
      Epithelial  = 1.8,
      Muscle      = 1.5
    ),
    alpha_subtypes = list(
      Immune      = 8,
      Stromal     = 8,
      Endothelial = 8,
      Epithelial  = 8,
      Muscle      = 8
    ),
    sigma_bio   = 0.15,
    sigma_tech  = 0.02,
    seed        = 1234
) {

  # ── 0. Load BlueCode from the package namespace ───────────────────────────
  # Avoids relying on the user having BlueCode in their global environment.
  W <- get("BlueCode", envir = asNamespace(utils::packageName()))

  # ── 1. Validate cell_groups against BlueCode column names ─────────────────
  cell_groups <- .bluecode_cell_groups
  stopifnot(
    "Cell types in .bluecode_cell_groups do not match colnames(BlueCode)" =
      setequal(unlist(cell_groups), colnames(W))
  )

  # ── 2. Reorder W columns to follow cell_groups compartment order ──────────
  ordered_celltypes <- unlist(cell_groups, use.names = FALSE)
  W <- W[, ordered_celltypes, drop = FALSE]

  # ── 3. Generate hierarchical proportions (real cell-type names) ────────────
  prop <- .generate_bluecode_prop(
    nSample        = nSample,
    cell_groups    = cell_groups,
    alpha_groups   = alpha_groups,
    alpha_subtypes = alpha_subtypes,
    seed           = seed
  )

  # ── 4. Linear mixing: B = W %*% t(prop) ───────────────────────────────────
  W_mat <- as.matrix(W)
  bulk  <- W_mat %*% t(prop)
  rownames(bulk) <- rownames(W_mat)
  colnames(bulk) <- rownames(prop)

  # ── 5. Biological noise: multiplicative log-normal ─────────────────────────
  # seed + 1L keeps proportion and noise sampling independently reproducible
  set.seed(seed + 1L)
  bio_noise  <- matrix(
    rnorm(length(bulk), mean = 0, sd = sigma_bio),
    nrow = nrow(bulk), ncol = ncol(bulk)
  )
  bulk_noisy <- bulk * exp(bio_noise)

  # ── 6. Technical noise: additive, scaled to global expression sd ───────────
  tech_noise <- matrix(
    rnorm(length(bulk), mean = 0, sd = sigma_tech * sd(bulk)),
    nrow = nrow(bulk), ncol = ncol(bulk)
  )
  bulk_noisy <- bulk_noisy + tech_noise

  # ── 7. Clip negative values ────────────────────────────────────────────────
  bulk_noisy[bulk_noisy < 0] <- 0

  rownames(bulk_noisy) <- rownames(bulk)
  colnames(bulk_noisy) <- colnames(bulk)

  # ── 8. Return standardised list ───────────────────────────────────────────
  return(list(
    p = as.data.frame(prop),
    W = as.data.frame(W),
    B = as.data.frame(bulk_noisy)
  ))
}
