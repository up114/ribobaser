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
