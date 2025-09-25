test_that("select_genes keeps dimensions and names consistent", {
  set.seed(1)
  G <- 6; S <- 4
  RIBO <- matrix(rpois(G*S, 50), nrow = G,
                 dimnames = list(paste0("g", 1:G), paste0("s", 1:S)))
  RNA  <- matrix(rpois(G*S, 60), nrow = G,
                 dimnames = list(paste0("g", 1:G), paste0("s", 1:S)))

  out <- select_genes(RIBO, RNA, cpm = 1, fraction = 0.5, dummy = FALSE)
  expect_true(is.matrix(out$ribo))
  expect_true(is.matrix(out$rna))
  expect_equal(colnames(out$ribo), colnames(RIBO))
  expect_equal(colnames(out$rna),  colnames(RNA))
  expect_equal(rownames(out$ribo), rownames(out$rna))
  expect_setequal(out$kept_genes, rownames(out$ribo))
})

test_that("select_genes removes low CPM genes and adds dummy row when requested", {
  # Construct data where g_low has tiny counts across all samples in both modalities
  RIBO <- rbind(
    g_keep1 = c(100, 120, 110),
    g_keep2 = c(200, 180, 190),
    g_low   = c(  1,   1,   1)
  )
  RNA <- rbind(
    g_keep1 = c(150, 160, 170),
    g_keep2 = c(300, 280, 260),
    g_low   = c(  0,   1,   0)
  )
  colnames(RIBO) <- colnames(RNA) <- paste0("s", 1:3)

  out <- select_genes(RIBO, RNA, cpm = 5, fraction = 0.8, dummy = TRUE)

  # g_low should be removed
  expect_false("g_low" %in% rownames(out$ribo))
  expect_false("g_low" %in% rownames(out$rna))
  expect_true("DUMMY_REMOVED" %in% rownames(out$ribo))
  expect_true("DUMMY_REMOVED" %in% rownames(out$rna))

  # Dummy counts equal the column sums of removed genes in each modality
  removed_ribo_colsums <- colSums(RIBO["g_low", , drop = FALSE])
  removed_rna_colsums  <- colSums(RNA ["g_low", , drop = FALSE])
  expect_equal(out$ribo["DUMMY_REMOVED", ], removed_ribo_colsums)
  expect_equal(out$rna ["DUMMY_REMOVED", ], removed_rna_colsums)

  # Kept gene set recorded correctly
  expect_setequal(out$kept_genes, c("g_keep1", "g_keep2"))
  expect_setequal(out$removed_genes, "g_low")
})

test_that("select_genes returns only dummy row if no genes pass", {
  # Everything will be below threshold
  RIBO <- matrix(1, nrow = 3, ncol = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:3)))
  RNA  <- matrix(1, nrow = 3, ncol = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:3)))

  out <- select_genes(RIBO, RNA, cpm = 1e6, fraction = 1, dummy = TRUE)
  expect_equal(rownames(out$ribo), "DUMMY_REMOVED")
  expect_equal(rownames(out$rna),  "DUMMY_REMOVED")

  # Dummy equals total col sums of all removed genes
  expect_equal(out$ribo["DUMMY_REMOVED", ], colSums(RIBO))
  expect_equal(out$rna ["DUMMY_REMOVED", ], colSums(RNA))
  expect_length(out$kept_genes, 0)
  expect_setequal(out$removed_genes, rownames(RIBO))
})

test_that("select_genes handles sample and gene intersections", {
  # Mismatched names but overlapping sets
  RIBO <- matrix(10, nrow = 3, ncol = 3,
                 dimnames = list(c("g1","g2","g3"), c("s1","s2","s3")))
  RNA  <- matrix(10, nrow = 4, ncol = 4,
                 dimnames = list(c("g0","g1","g3","g9"), c("s0","s1","s2","s3")))

  out <- select_genes(RIBO, RNA, cpm = 1, fraction = 0.5, dummy = FALSE)

  # Only shared genes and samples remain
  expect_setequal(rownames(out$ribo), c("g1","g3"))
  expect_setequal(colnames(out$ribo), c("s1","s2","s3"))
  expect_equal(rownames(out$ribo), rownames(out$rna))
  expect_equal(colnames(out$ribo), colnames(out$rna))
})

test_that("select_genes errors cleanly when no shared genes or samples", {
  R1 <- matrix(10, nrow = 2, ncol = 2,
               dimnames = list(c("a","b"), c("s1","s2")))
  R2 <- matrix(10, nrow = 2, ncol = 2,
               dimnames = list(c("x","y"), c("t1","t2")))

  expect_error(select_genes(R1, R2), "No shared genes")
  # tweak so genes match but samples do not
  rownames(R2) <- rownames(R1)
  expect_error(select_genes(R1, R2), "No shared samples")
})

test_that("select_genes is stable with single kept gene due to drop = FALSE", {
  RIBO <- rbind(
    g_keep = c(100, 100, 100),
    g_low  = c(  1,   1,   1)
  )
  RNA <- rbind(
    g_keep = c(200, 200, 200),
    g_low  = c(  0,   1,   0)
  )
  colnames(RIBO) <- colnames(RNA) <- paste0("s", 1:3)

  out <- select_genes(RIBO, RNA, cpm = 10, fraction = 0.9, dummy = FALSE)
  expect_true(is.matrix(out$ribo))
  expect_true(is.matrix(out$rna))
  expect_equal(rownames(out$ribo), "g_keep")
  expect_equal(dim(out$ribo), c(1, 3))
})
