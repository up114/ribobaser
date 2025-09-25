# Tests for zero_imputation

test_that("zero_imputation replaces zeros with ones", {
  RIBO <- matrix(
    c(0, 2, 0,
      4, 5, 6),
    nrow = 2, byrow = T,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )
  RNA <- matrix(
    c(7, 8, 9,
      0, 0, 3),
    nrow = 2, byrow = T,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )

  res <- zero_imputation(RIBO, RNA)

  expect_true(all(res$ribo >= 1))
  expect_true(all(res$rna  >= 1))
  expect_equal(res$ribo[1, 1], 1)
  expect_equal(res$ribo[1, 3], 1)
  expect_equal(res$rna [2, 1], 1)
  expect_equal(res$rna [2, 2], 1)
})

test_that("zero_imputation preserves non-zero values and structure", {
  RIBO <- data.frame(
    s1 = c(0, 5),
    s2 = c(2, 0),
    row.names = c("g1", "g2")
  )
  RNA <- data.frame(
    s1 = c(3, 4),
    s2 = c(0, 6),
    row.names = c("g1", "g2")
  )

  res <- zero_imputation(RIBO, RNA)

  expect_identical(dim(res$ribo), dim(as.matrix(RIBO)))
  expect_identical(dim(res$rna),  dim(as.matrix(RNA)))
  expect_equal(res$ribo["g2", "s1"], 5)
  expect_equal(res$rna ["g2", "s2"], 6)
})
