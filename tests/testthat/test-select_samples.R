# Tests for select_samples

test_that("select_samples returns inputs unchanged when filtering disabled", {
  RIBO <- matrix(
    c(10, 20, 30,
      40, 50, 60),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )
  RNA <- RIBO + 5

  res <- select_samples(RIBO, RNA, high_dummy_percentage = 0, min_r2 = 0)

  expect_identical(res$ribo, RIBO)
  expect_identical(res$rna, RNA)
})

test_that("select_samples validates high_dummy_percentage input", {
  RIBO <- matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))
  RNA  <- RIBO

  expect_error(
    select_samples(RIBO, RNA, high_dummy_percentage = -0.1),
    "high_dummy_percentage must be a number between 0 and 1"
  )
  expect_error(
    select_samples(RIBO, RNA, high_dummy_percentage = c(0.2, 0.3)),
    "high_dummy_percentage must be a number between 0 and 1"
  )
})

test_that("select_samples validates min_r2 input", {
  RIBO <- matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))
  RNA  <- RIBO

  expect_error(
    select_samples(RIBO, RNA, min_r2 = -0.1),
    "min_r2 must be a number between 0 and 1"
  )
  expect_error(
    select_samples(RIBO, RNA, min_r2 = c(0.1, 0.2)),
    "min_r2 must be a number between 0 and 1"
  )
})

test_that("select_samples validates min_periodicity input", {
  RIBO <- matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))
  RNA  <- RIBO

  expect_error(
    select_samples(RIBO, RNA, min_periodicity = -0.1),
    "min_periodicity must be a number between 0 and 1"
  )
  expect_error(
    select_samples(RIBO, RNA, min_periodicity = c(0.1, 0.2)),
    "min_periodicity must be a number between 0 and 1"
  )
})

test_that("select_samples errors when dummy row missing but filtering requested", {
  RIBO <- matrix(
    c(100, 150,
      200, 250),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  RNA <- RIBO

  expect_error(
    select_samples(RIBO, RNA, high_dummy_percentage = 0.2),
    "DUMMY_REMOVED row must exist"
  )
})

test_that("select_samples errors when periodicity metadata missing", {
  RIBO <- matrix(1, nrow = 2, ncol = 1, dimnames = list(paste0("g", 1:2), "unknown_sample"))
  RNA  <- RIBO

  expect_error(
    select_samples(RIBO, RNA, min_periodicity = 0.5),
    "Missing periodicity scores"
  )
})

test_that("select_samples removes samples exceeding dummy threshold in either modality", {
  RIBO <- rbind(
    g1 = c(100,  50, 100,  70, 100),
    g2 = c(200,  30,  90,  50, 200),
    DUMMY_REMOVED = c(10, 100, 10, 45, 90)
  )
  RNA <- rbind(
    g1 = c(150, 60,  80, 60, 100),
    g2 = c(160, 40,  50, 40, 200),
    DUMMY_REMOVED = c(5, 110, 100, 15, 90)
  )
  colnames(RIBO) <- colnames(RNA) <- paste0("s", 1:5)

  res <- select_samples(RIBO, RNA, high_dummy_percentage = 0.3)

  expect_identical(colnames(res$ribo), c("s1", "s5"))
  expect_identical(colnames(res$rna),  c("s1", "s5"))
  expect_equal(res$ribo[, "s5"], RIBO[, "s5"])
  expect_equal(res$rna[,  "s5"], RNA[,  "s5"])
})

test_that("select_samples removes samples with low R-squared correspondence", {
  skip_if_not_installed("compositions")

  RIBO <- rbind(
    g1 = c(100,  50,  70),
    g2 = c(200,  40,  20),
    g3 = c(300,  60, 500)
  )
  RNA <- rbind(
    g1 = c(110, 300,  3),
    g2 = c(210, 100,  25),
    g3 = c(290, 250, 800)
  )
  colnames(RIBO) <- colnames(RNA) <- paste0("s", 1:3)

  # Sample s1 has strong correspondence, s2 low, s3 moderate
  res <- select_samples(RIBO, RNA, min_r2 = 0.7)

  expect_identical(colnames(res$ribo), "s1")
  expect_identical(colnames(res$rna),  "s1")
})

test_that("select_samples removes samples below periodicity threshold", {
  data("Ribobase_QC_non_dedup_data")
  low_idx <- which(Ribobase_QC_non_dedup_data$`Periodicity score` < 0.5)[1:2]
  high_idx <- which(Ribobase_QC_non_dedup_data$`Periodicity score` > 0.7)[1:2]
  exps <- Ribobase_QC_non_dedup_data$Experiment[c(low_idx, high_idx)]

  RIBO <- matrix(
    c(10, 20,  30,  40,
      50, 60,  70,  80,
      90, 100, 110, 120),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("g", 1:3), exps)
  )
  RNA <- RIBO + 5

  res <- select_samples(RIBO, RNA, min_periodicity = 0.6)

  expect_false(any(colnames(res$ribo) %in% exps[1:2]))
  expect_setequal(colnames(res$ribo), exps[3:4])
  expect_setequal(colnames(res$rna),  exps[3:4])
})

test_that("select_samples applies both filters", {
  skip_if_not_installed("compositions")

  RIBO <- rbind(
    g1 = c(100,  50, 400),
    g2 = c(150,  40, 300),
    DUMMY_REMOVED = c(10, 150, 20)
  )
  RNA <- rbind(
    g1 = c(105, 150, 5),
    g2 = c(155, 130, 320),
    DUMMY_REMOVED = c(5, 200, 15)
  )
  colnames(RIBO) <- colnames(RNA) <- paste0("s", 1:3)

  res <- select_samples(RIBO, RNA, high_dummy_percentage = 0.3, min_r2 = 0.8)
  res
  # s2 removed by dummy filter, s3 removed by low R^2
  expect_identical(colnames(res$ribo), "s1")
  expect_identical(colnames(res$rna),  "s1")
})

test_that("select_samples returns empty matrices when all samples removed", {
  skip_if_not_installed("compositions")

  RIBO <- rbind(
    g1 = c(10,  20),
    DUMMY_REMOVED = c(5,  6)
  )
  RNA <- rbind(
    g1 = c(12, 18),
    DUMMY_REMOVED = c(6,  7)
  )
  colnames(RIBO) <- colnames(RNA) <- c("s1", "s2")

  res <- select_samples(RIBO, RNA, high_dummy_percentage = 0.2, min_r2 = 0.9)

  expect_true(is.matrix(res$ribo))
  expect_true(is.matrix(res$rna))
  expect_identical(rownames(res$ribo), rownames(RIBO))
  expect_identical(rownames(res$rna),  rownames(RNA))
  expect_equal(ncol(res$ribo), 0)
  expect_equal(ncol(res$rna),  0)
})
