#' @title BlueCode Reference Signature Matrix
#'
#' @description Gene expression reference matrix for bulk RNA-seq deconvolution,
#' derived from sorted cell populations spanning Immune, Stromal, Endothelial,
#' epithelial, and muscle compartments. BlueCode serves as the default
#' reference for dicepro and directly informs the hierarchical Dirichlet
#' simulation framework.
#'
#' @docType data
#' @keywords data
#' @name BlueCode
#' @usage data(BlueCode)
#'
#' @format A numeric matrix of dimensions \eqn{G \times 34}, where \eqn{G}
#' denotes the number of genes after quality filtering. Rows are labeled
#' by HGNC gene symbol; columns are labeled by cell-type name. All values
#' are non-negative read counts prior to normalization; apply
#' \code{.normalize_zscore_per_gene} or TPM scaling before use in
#' downstream analyses.
#'
#' @details
#' ## Construction
#'
#' BlueCode was built from publicly available sorted bulk RNA-seq profiles
#' downloaded from the ENCODE and GEO repositories. For each cell population,
#' replicate profiles were averaged after quantile normalization. Genes with
#' zero variance across all 34 profiles, or with missing values in more than
#' 10 \% of replicates, were removed. The resulting matrix was \emph{not}
#' log-transformed so that it can be used directly with both Gaussian and
#' Poisson deconvolution models.
#'
#' ## Compartment structure
#'
#' The 34 cell types are partitioned into five major tissue compartments.
#' The column ordering in the matrix follows this partition exactly (columns
#' 1--9, 10--17, 18--20, 21--25, 26--34).
#'
#' \strong{Immune compartment — columns 1--9:}
#' \enumerate{
#'   \item \code{B.cell.naive}
#'   \item \code{B.cell.memory}
#'   \item \code{T.Cell.CD4}
#'   \item \code{T.cell.CD4.memory}
#'   \item \code{T.cell.CD8.memory}
#'   \item \code{NK.cell.CD56}
#'   \item \code{Monocyte.CD14}
#'   \item \code{Macrophage}
#'   \item \code{Mature.neutrophil}
#' }
#'
#' \strong{Stromal compartment — columns 10--17:}
#' \enumerate{
#'   \item \code{Fibroblast.cardiac.ventricle}
#'   \item \code{Fibroblast.of.arm}
#'   \item \code{Fibroblast.of.lung}
#'   \item \code{Fibroblast.of.the.aortic.adventitia}
#'   \item \code{Fibroblast.or.papilla.dermal.cell}
#'   \item \code{MSC.like.pluripotent.cell}
#'   \item \code{Chondrocyte.articular}
#'   \item \code{Osteoblast}
#' }
#'
#' \strong{Endothelial compartment — columns 18--20:}
#' \enumerate{
#'   \item \code{Endothelial.large.blood.vessel}
#'   \item \code{Endothelial.microvascular.mammary.or.endometrial}
#'   \item \code{Endothelial.microvascular.non.reproductive}
#' }
#'
#' \strong{Epithelial compartment — columns 21--25:}
#' \enumerate{
#'   \item \code{Epithelial.cell.mammary}
#'   \item \code{Epithelial.renal.cortical}
#'   \item \code{Epithelial.respiratory}
#'   \item \code{Hair.follicular.keratinocyte}
#'   \item \code{Melanocyte.of.skin}
#' }
#'
#' \strong{Muscle compartment — columns 26--34:}
#' \enumerate{
#'   \item \code{Smooth.muscle.cell.aortic}
#'   \item \code{Smooth.muscle.cell.bronchial}
#'   \item \code{Smooth.muscle.cell.coronary.artery}
#'   \item \code{Smooth.muscle.cell.pulmonary.artery}
#'   \item \code{Smooth.muscle.cell.of.bladder}
#'   \item \code{Smooth.muscle.cell.of.trachea}
#'   \item \code{Smooth.muscle.cell.uterine}
#'   \item \code{Myocyte.regular.cardiac}
#'   \item \code{Myometrial.cell}
#' }
#'
#' ## Link to the simulation framework
#'
#' This compartment structure directly informs the hierarchical Dirichlet
#' model in \code{generateProp} with \code{scenario = "hierarchical"}.
#' The Dirichlet concentration parameters at each level
#' (\eqn{\alpha_\text{Immune} = 6}, \eqn{\alpha_\text{Stromal} = 2},
#' \eqn{\alpha_\text{Endothelial} = \alpha_\text{Epithelial} = 1.5},
#' \eqn{\alpha_\text{Muscle} = 1}) reflect the expected relative abundance
#' of each compartment in mixed tissue samples, so that simulated proportions
#' are biologically realistic with respect to BlueCode cell-type coverage.
#'
#' ## Recommended usage
#'
#' ```r
#' data(BlueCode)
#'
#' # Quick inspection
#' dim(BlueCode)          # G x 34
#' colnames(BlueCode)     # 34 cell-type labels
#'
#' # Use directly in dicepro
#' out <- dicepro(
#'   reference    = BlueCode,
#'   bulk         = my_bulk_matrix,
#'   methodDeconv = "FARDEEP"
#' )
#' ```
#'
#' @source
#' BlueCode reference matrix constructed for the dicepro benchmarking
#' framework from ENCODE (\url{https://www.encodeproject.org}) and GEO
#' (\url{https://www.ncbi.nlm.nih.gov/geo/}) sorted RNA-seq profiles.
#' See the accompanying manuscript for full details on data sources,
#' replicate averaging, quality filtering, and normalization.
"BlueCode"


#' @title CellMixtures Bulk RNA-seq Dataset
#'
#' @description Experimentally constructed bulk RNA-seq mixtures of sorted
#' cell populations, designed for benchmarking deconvolution algorithms.
#' Each sample is a known mixture of cell types drawn from the BlueCode
#' reference panel, making ground-truth proportions available for
#' quantitative evaluation.
#'
#' @docType data
#' @keywords data
#' @name CellMixtures
#' @usage data(CellMixtures)
#'
#' @format A numeric matrix of dimensions \eqn{G \times 12}, where \eqn{G}
#' denotes the number of genes (~31 400 after quality filtering). Rows are
#' labeled by HGNC gene symbol; columns are labeled \code{A} through
#' \code{L}, corresponding to 12 experimentally mixed samples. Values are
#' raw read counts prior to normalization.
#'
#' @details
#' ## Experimental design
#'
#' The 12 mixture samples (A--L) were prepared by combining sorted cell
#' populations at known proportions across the five tissue compartments
#' represented in \code{BlueCode}. Each sample targets a distinct
#' mixture composition, providing a controlled benchmark spanning a wide
#' range of cell-type abundance profiles — from immune-dominated samples
#' to Stromal- or muscle-enriched configurations.
#'
#' ## Gene coverage
#'
#' CellMixtures covers approximately 31 400 gene symbols, a super-set of
#' those present in BlueCode. dicepro automatically intersects the two
#' gene sets before deconvolution; the effective number of informative
#' genes is therefore determined by the BlueCode reference (see
#' \code{BlueCode} for details on its gene filtering criteria).
#' A gene-overlap check prior to running \code{dicepro} is
#' recommended:
#'
#' ```r
#' data(BlueCode)
#' data(CellMixtures)
#' n_common <- length(intersect(rownames(BlueCode), rownames(CellMixtures)))
#' cat(sprintf("Common genes: %d / %d reference genes\n",
#'             n_common, nrow(BlueCode)))
#' ```
#'
#' ## normalization
#'
#' Raw counts should be normalized before deconvolution.
#' \code{dicepro} applies \code{.normalize_zscore_per_gene} internally;
#' no preprocessing is required when using the main function.
#' For exploratory analyses, log-CPM or TPM normalization is recommended:
#'
#' ```r
#' # log2-CPM normalization
#' cpm <- sweep(CellMixtures, 2, colSums(CellMixtures) / 1e6, FUN = "/")
#' log2_cpm <- log2(cpm + 1)
#' ```
#'
#' ## Recommended usage
#'
#' ```r
#' data(BlueCode)
#' data(CellMixtures)
#'
#' out <- dicepro(
#'   reference             = BlueCode,
#'   bulk                  = CellMixtures,
#'   methodDeconv          = "FARDEEP",
#'   bulkName              = "CellMixtures",
#'   refName               = "BlueCode",
#'   hp_max_evals          = 100L,
#'   hspaceTechniqueChoose = "all"
#' )
#'
#' # Inspect estimated proportions
#' head(out$H)
#' ```
#'
#' ## Relationship to BlueCode
#'
#' CellMixtures and BlueCode are designed as a paired benchmark: the
#' reference matrix is derived from the same cell populations used to
#' construct the mixtures, ensuring that all cell types present in the bulk
#' are represented in the reference. This controlled setting enables direct
#' evaluation of the reconstruction accuracy of dicepro.
#'
#' @source
#' CellMixtures was constructed for the DICEPro benchmarking framework.
#' Sorted cell populations were obtained from ENCODE
#' (\url{https://www.encodeproject.org}) and mixed \emph{in silico} at
#' predefined proportions from experimentally validated RNA-seq profiles.
#' See the accompanying manuscript for full details on mixture design,
#' library preparation, and sequencing.
"CellMixtures"
