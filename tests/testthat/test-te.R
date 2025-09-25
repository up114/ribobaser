# Tests for translation efficiency estimator
test_that("te returns near-zero residuals when RIBO and RNA match", {
  skip_if_not_installed("propr")
  skip_if_not_installed("compositions")
  skip_if_not_installed("foreach")

  set.seed(123)
  samples <- paste0("s", 1:3)
  genes <- paste0("g", 1:4)
  counts <- matrix(
    runif(length(samples) * length(genes), min = 50, max = 200),
    nrow = length(samples),
    dimnames = list(samples, genes)
  )

  out <- te(counts, counts)

  expect_true(is.matrix(out))
  expect_equal(rownames(out), genes)
  expect_equal(colnames(out), samples)

  zero_matrix <- matrix(
    0,
    nrow = length(genes),
    ncol = length(samples),
    dimnames = list(genes, samples)
  )
  expect_equal(out, zero_matrix, tolerance = 1e-6)
})

test_that("te aligns overlapping samples and genes", {
  skip_if_not_installed("propr")
  skip_if_not_installed("compositions")
  skip_if_not_installed("foreach")

  RIBO <- matrix(
    c(
      100, 120, 140, 160,
      110, 130, 150, 170,
      120, 140, 160, 180
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("s1", "s2", "s3"),
      c("g1", "g2", "g3", "g4")
    )
  )
  RNA <- matrix(
    c(
      90, 110, 130, 150,
      95, 115, 135, 155,
      105, 125, 145, 165
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("s2", "s3", "s4"),
      c("g0", "g2", "g3", "g5")
    )
  )

  out <- te(RIBO, RNA)

  expect_true(is.matrix(out))
  expect_setequal(rownames(out), c("g2", "g3"))
  expect_setequal(colnames(out), c("s2", "s3"))
  expect_equal(dim(out), c(2, 2))
})

test_that("te errors when there are no shared samples or genes", {
  skip_if_not_installed("propr")
  skip_if_not_installed("compositions")
  skip_if_not_installed("foreach")

  R1 <- matrix(
    c(
      10, 20,
      30, 40
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("s1", "s2"),
      c("g1", "g2")
    )
  )
  R2 <- matrix(
    c(
      5, 15,
      25, 35
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("x1", "x2"),
      c("g1", "g2")
    )
  )

  expect_error(te(R1, R2), "No shared samples")

  rownames(R2) <- rownames(R1)
  colnames(R2) <- c("g9", "g10")

  expect_error(te(R1, R2), "No shared genes")
})

test_that("te logratio method matches manual clr difference", {
  RIBO <- matrix(
    c(
      20, 30, 40,
      35, 45, 55
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("s1", "s2"),
      c("g1", "g2", "g3")
    )
  )
  RNA <- matrix(
    c(
      10, 20, 30,
      25, 35, 45
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("s1", "s2"),
      c("g1", "g2", "g3")
    )
  )

  clr_manual <- function(mat) {
    log_mat <- log(mat)
    sweep(log_mat, 1, rowMeans(log_mat), "-")
  }

  expected <- t(clr_manual(RIBO) - clr_manual(RNA))
  out <- te(RIBO, RNA, method = "logratio")

  expect_true(is.matrix(out))
  expect_equal(rownames(out), c("g1", "g2", "g3"))
  expect_equal(colnames(out), c("s1", "s2"))
  expect_equal(out, expected)
})

test_that("te logratio errors on non positive inputs", {
  RIBO <- matrix(
    c(
      10, 0,
      15, 25
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("s1", "s2"),
      c("g1", "g2")
    )
  )
  RNA <- matrix(
    c(
      12, 16,
      18, 20
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("s1", "s2"),
      c("g1", "g2")
    )
  )

  expect_error(te(RIBO, RNA, method = "logratio"), "strictly positive")
})
