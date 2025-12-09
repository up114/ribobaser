# Tests for aggregate_samples

test_that("aggregate_samples computes study means", {
  data("Ribobase_QC_dedup_data")
  studies <- Ribobase_QC_dedup_data$Experiment[1:4]

  TE <- matrix(
    c(5, 7, 3, 9,
      2, 4, 6, 8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), studies)
  )

  agg <- aggregate_samples(TE, group_by = "study")

  expect_true(is.matrix(agg))
  expect_identical(rownames(agg), rownames(TE))
  expect_identical(colnames(agg), unique(Ribobase_QC_dedup_data$Study[match(studies, Ribobase_QC_dedup_data$Experiment)]))

  expected_first <- rowMeans(TE[, 1:3, drop = FALSE])
  expected_second <- TE[, 4]
  expect_equal(agg[, 1], expected_first)
  expect_equal(agg[, 2], expected_second)
})

test_that("aggregate_samples computes cell line means", {
  data("Ribobase_QC_dedup_data")
  exps <- c(Ribobase_QC_dedup_data$Experiment[1:3], Ribobase_QC_dedup_data$Experiment[6])

  TE <- matrix(
    c(10, 12, 14, 20,
      5,   7,  9, 25),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), exps)
  )

  agg <- aggregate_samples(TE, group_by = "cell_line")

  expect_identical(colnames(agg), c("HEK293", "RPE-1"))

  expected_hek <- rowMeans(TE[, 1:3, drop = FALSE])
  expected_rpe <- TE[, 4]
  expect_equal(agg[, "HEK293"], expected_hek)
  expect_equal(agg[, "RPE-1"], expected_rpe)
})

test_that("aggregate_samples errors on unknown samples", {
  TE <- matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "unknown_sample"))
  expect_error(aggregate_samples(TE), "Missing metadata")
})

test_that("aggregate_samples gives same result for string and data.frame metadata (study)", {
  data("Ribobase_QC_dedup_data")

  studies <- Ribobase_QC_dedup_data$Experiment[1:4]

  TE <- matrix(
    c(5, 7, 3, 9,
      2, 4, 6, 8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), studies)
  )

  agg_str <- aggregate_samples(
    TE,
    metadata = "Ribobase_QC_dedup_data",
    group_by = "study",
    fun = "mean"
  )

  agg_df <- aggregate_samples(
    TE,
    metadata = Ribobase_QC_dedup_data,
    group_by = "study",
    fun = "mean"
  )

  expect_identical(agg_str, agg_df)
})

test_that("aggregate_samples works with a custom metadata data.frame", {
  # tiny fake metadata
  metadata_custom <- data.frame(
    Experiment = c("S1", "S2", "S3", "S4"),
    Study      = c("GSE_A", "GSE_A", "GSE_B", "GSE_B"),
    `Cell line` = c("CL1", "CL1", "CL2", "CL2"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  TE <- matrix(
    c(1, 3, 5, 7,
      2, 4, 6, 8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), metadata_custom$Experiment)
  )

  agg <- aggregate_samples(
    TE,
    metadata = metadata_custom,
    group_by = "study",
    fun = "mean"
  )

  # Expected: two studies: GSE_A (S1,S2) and GSE_B (S3,S4)
  expect_identical(colnames(agg), c("GSE_A", "GSE_B"))

  expected_A <- rowMeans(TE[, c("S1", "S2"), drop = FALSE])
  expected_B <- rowMeans(TE[, c("S3", "S4"), drop = FALSE])

  expect_equal(agg[, "GSE_A"], expected_A)
  expect_equal(agg[, "GSE_B"], expected_B)
})

