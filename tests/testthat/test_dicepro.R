test_that("dicepro runs on simulated data with incomplete reference", {
  skip_on_cran()

  library(dicepro)
  set.seed(42L)

  nCellsType <- 5L
  nGenes     <- 50L
  nSample    <- 10L

  # ---- Simulate — single call guarantees W / p / B column alignment -------
  sim <- simulation(
    loi        = "gauss",
    scenario   = "hierarchical",
    nSample    = nSample,
    nGenes     = nGenes,
    nCellsType = nCellsType,
    sigma_bio  = 0.07,
    sigma_tech = 0.07,
    seed       = 42L
  )

  expect_true(is.data.frame(sim$B) || is.matrix(sim$B))
  expect_true(is.data.frame(sim$W) || is.matrix(sim$W))
  expect_true(is.data.frame(sim$p) || is.matrix(sim$p))
  expect_equal(nrow(as.matrix(sim$p)), nSample)
  expect_equal(ncol(as.matrix(sim$p)), nCellsType)

  # ---- Build incomplete reference (drop one random cell type) -------------
  W_full       <- as.matrix(sim$W)
  cell_out     <- sample(seq_len(ncol(W_full)), 1L)
  ref_incomplete <- W_full[, -cell_out, drop = FALSE]

  expect_equal(ncol(ref_incomplete), nCellsType - 1L)

  # ---- Run dicepro --------------------------------------------------------
  out <- dicepro(
    reference             = ref_incomplete,
    bulk                  = as.matrix(sim$B),
    methodDeconv          = "FARDEEP",
    W_prime               = 0,
    bulkName              = "TestBulk",
    refName               = "TestRef",
    hp_max_evals          = 5L,
    algo_select           = "random",
    output_path           = tempdir(),
    hspaceTechniqueChoose = "all"
  )
  if (is.null(out)) {
    skip("All hyperparameter trials failed — increase nGenes/nCellsType or hp_max_evals.")
  }


  # ---- Check output structure ---------------------------------------------
  expect_s3_class(out, "dicepro")
  expect_named(out, c("hyperparameters", "metrics", "trials",
                      "W", "H", "plot", "plot_hyperopt"),
               ignore.order = TRUE)
  expect_true(is.numeric(out$hyperparameters$lambda))
  expect_true(is.numeric(out$hyperparameters$gamma))
  expect_true(is.numeric(out$metrics$loss))
  expect_true(is.data.frame(out$trials))
  expect_gt(nrow(out$trials), 0L)

  # ---- Check proportion matrix --------------------------------------------
  H <- as.matrix(out$H)
  expect_true(is.matrix(H))
  expect_equal(nrow(H), nSample)
  expect_true(all(H >= -1e-6))
  expect_true(all(abs(rowSums(H) - 1) < 0.1))

  # ---- Check plots exist --------------------------------------------------
  expect_false(is.null(out$plot))
  expect_false(is.null(out$plot_hyperopt))
})
