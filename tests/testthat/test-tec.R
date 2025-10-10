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

  expect_equal(rho_mat, expected)
})

test_that("tec errors on invalid inputs", {
  skip_if_not_installed("propr")

  bad_te <- matrix(c(NA_real_, 0.1, 0.2, 0.3), nrow = 2)
  expect_error(tec(bad_te), "without missing data")

  too_small <- matrix(0.5, nrow = 1, ncol = 3)
  expect_error(tec(too_small), "at least two genes")
})
