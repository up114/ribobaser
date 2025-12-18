test_that("tec matches propr:::lr2rho for regression TE", {
  skip_if_not_installed("propr")
  skip_if_not_installed("compositions")
  skip_if_not_installed("foreach")

  data("ribo_raw_human_cap_995", package = "ribobaser")
  data("rnaseq_raw_human_cap_995", package = "ribobaser")

  genes <- rownames(ribo_raw_human_cap_995)[1:4]
  samples <- colnames(ribo_raw_human_cap_995)[1:6]

  ribo <- t(as.matrix(ribo_raw_human_cap_995[genes, samples])) + 1
  rna  <- t(as.matrix(rnaseq_raw_human_cap_995[genes, samples])) + 1

  te_mat <- te(ribo, rna, method = "regression")
  rho_mat <- tec(te_mat)

  expect_true(is.matrix(rho_mat))
  expect_equal(rownames(rho_mat), genes)
  expect_equal(colnames(rho_mat), genes)

  rho_direct <- propr:::lr2rho(t(te_mat))
  if (is.matrix(rho_direct)) {
    expected <- rho_direct
  } else {
    pair_count <- length(rho_direct)
    expected <- matrix(1, nrow = length(genes), ncol = length(genes))
    expected[upper.tri(expected)] <- rho_direct
    expected[lower.tri(expected)] <- t(expected)[lower.tri(expected)]
  }
  dimnames(expected) <- list(genes, genes)

  expect_equal(rho_mat, expected)
})

test_that("tec errors on invalid inputs", {
  skip_if_not_installed("propr")

  bad_te <- matrix(c(NA_real_, 0.1, 0.2, 0.3), nrow = 2)
  expect_error(tec(bad_te), "without missing data")

  too_small <- matrix(0.5, nrow = 1, ncol = 3)
  expect_error(tec(too_small), "at least two genes")
})

test_that("tec(wgcna) returns symmetric TOM with expected dimnames", {
  skip_if_not_installed("WGCNA")

  set.seed(123)
  genes   <- paste0("gene", 1:6)
  samples <- paste0("sample", 1:10)

  # genes x samples with non-zero variance and no NAs
  TE <- matrix(
    rnorm(length(genes) * length(samples)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )

  tom <- tec(TE, method = "wgcna")

  expect_true(is.matrix(tom))
  expect_equal(dim(tom), c(length(genes), length(genes)))
  expect_equal(rownames(tom), genes)
  expect_equal(colnames(tom), genes)

  # TOM should be symmetric
  expect_equal(tom, t(tom), tolerance = 1e-8)
})

test_that("tec(glasso) returns symmetric precision matrix with expected dimnames", {
  skip_if_not_installed("huge")

  set.seed(123)
  genes   <- paste0("gene", 1:5)
  samples <- paste0("sample", 1:12)

  # genes x samples, random values ensure non-zero variance
  TE <- matrix(
    rnorm(length(genes) * length(samples)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )

  theta <- tec(TE, method = "glasso")

  expect_true(is.matrix(theta))
  expect_equal(dim(theta), c(length(genes), length(genes)))
  expect_equal(rownames(theta), genes)
  expect_equal(colnames(theta), genes)

  # precision matrix should be symmetric after your symmetrization
  expect_equal(theta, t(theta), tolerance = 1e-8)

  # diagonal typically positive (weak sanity check, not mathematically strict)
  expect_true(all(diag(theta) > 0))
})

test_that("tec(genie3) returns symmetric, reproducible weight matrix", {
  skip_on_cran()
  skip_if_not_installed("GENIE3")

  genes   <- paste0("gene", 1:5)
  samples <- paste0("sample", 1:8)

  set.seed(42)
  TE <- matrix(
    rnorm(length(genes) * length(samples)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )

  W1 <- tec(TE, method = "genie3")
  W2 <- tec(TE, method = "genie3")

  expect_true(is.matrix(W1))
  expect_equal(dim(W1), c(length(genes), length(genes)))
  expect_equal(rownames(W1), genes)
  expect_equal(colnames(W1), genes)

  # symmetric after (W + t(W))/2
  expect_equal(W1, t(W1), tolerance = 1e-8)

  # fixed seed inside tec() should make results reproducible
  expect_equal(W1, W2)
})

test_that("tec(glasso) drops near-zero variance genes", {
  skip_if_not_installed("huge")

  genes   <- c("stable1", "stable2", "variable1", "variable2")
  samples <- paste0("sample", 1:10)

  # two almost-constant genes, two variable genes
  stable_block <- matrix(1, nrow = 2, ncol = length(samples))
  variable_block <- matrix(
    rnorm(2 * length(samples)),
    nrow = 2
  )

  TE <- rbind(stable_block, variable_block)
  rownames(TE) <- genes
  colnames(TE) <- samples

  Theta <- tec(TE, method = "glasso")

  # Only variable genes should remain
  expect_true(all(rownames(Theta) %in% c("variable1", "variable2")))
  expect_true(all(colnames(Theta) %in% c("variable1", "variable2")))
  expect_equal(dim(Theta), c(2L, 2L))
})
