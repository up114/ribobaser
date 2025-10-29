# tests/testthat/test-pc_remove.R

test_that("pc_remove: rejects invalid inputs and missing names", {
  set.seed(1)
  M <- matrix(rnorm(20), nrow = 5, ncol = 4)                # no dimnames
  expect_error(pc_remove(M, method = "global"),
               "rownames.*colnames", ignore.case = TRUE)

  rownames(M) <- paste0("g", 1:5)
  colnames(M) <- paste0("s", 1:4)

  # non-matrix / non-numeric
  expect_error(pc_remove(as.data.frame(M), method = "global"))

  # NA/Inf check
  M_bad <- M
  M_bad[1, 1] <- NA_real_
  expect_error(pc_remove(M_bad, method = "global"),
               "NA/Inf", ignore.case = TRUE)
})

test_that("pc_remove: global k=0 performs per-sample demeaning", {
  set.seed(2)
  G <- 10; S <- 6
  M <- matrix(rnorm(G * S), nrow = G, ncol = S)
  rownames(M) <- paste0("g", seq_len(G))
  colnames(M) <- paste0("s", seq_len(S))

  out <- pc_remove(M, method = "global", n_pcs = 0L)
  X   <- out$te_corrected

  # When k=0, each sample (column) is demeaned across genes:
  # colMeans across *all genes kept for SVD* should be ~0.
  # Since we provided random data (non-constant rows), all rows are "kept".
  expect_equal(as.numeric(colMeans(X)), rep(0, S), tolerance = 1e-10)

  # dims and names preserved
  expect_identical(dim(X), dim(M))
  expect_identical(rownames(X), rownames(M))
  expect_identical(colnames(X), colnames(M))
  # n_pcs should be 0 (global integer)
  expect_identical(out$n_pcs, 0L)
})

test_that("pc_remove: zero-variance rows are unchanged and reported", {
  set.seed(3)
  G <- 8; S <- 6
  M <- matrix(rnorm(G * S), nrow = G, ncol = S)
  rownames(M) <- paste0("g", seq_len(G))
  colnames(M) <- paste0("s", seq_len(S))

  # Make one constant row
  M[3, ] <- 5

  out <- pc_remove(M, method = "global", n_pcs = 1L)
  X   <- out$te_corrected

  # The constant row must be unchanged
  expect_equal(X[3, ], M[3, ], tolerance = 0)

  # Row should appear in dropped_genes
  expect_true("g3" %in% out$dropped_genes)
})

test_that("pc_remove: k is clamped to valid range", {
  set.seed(4)
  G <- 7; S <- 5
  M <- matrix(rnorm(G * S), nrow = G, ncol = S)
  rownames(M) <- paste0("g", seq_len(G))
  colnames(M) <- paste0("s", seq_len(S))

  # Asking for an absurdly large k should clamp to min(n-1, p-1) in sample space:
  # Here Y is S x G, so max PCs = min(S, G) - 1 = min(5,7)-1 = 4 - actually 4 (because min-1).
  out <- pc_remove(M, method = "global", n_pcs = 1000L)
  expect_equal(out$n_pcs, min(S, G) - 1L)
})

test_that("pc_remove: removes a rank-1 sample-level signal when k=1", {
  set.seed(5)
  # Build Y (samples x genes) with a rank-1 signal + noise
  S <- 12; G <- 30
  s <- scale(rnorm(S), center = TRUE, scale = FALSE)  # sample signal (demeaned)
  g <- rnorm(G)                                       # gene weights
  Y <- s %*% t(g) + matrix(rnorm(S * G, sd = 0.1), nrow = S)

  # te_mat is genes x samples:
  M <- t(Y)
  rownames(M) <- paste0("g", seq_len(G))
  colnames(M) <- paste0("s", seq_len(S))

  # Remove 1 PC
  out <- pc_remove(M, method = "global", n_pcs = 1L)
  X   <- out$te_corrected  # genes x samples
  Yadj <- t(X)             # samples x genes

  # After removing the top sample PC, projection of Yadj onto the original s should be ~0
  # Compute correlation between s and each gene's profile; expect near 0 on average.
  cors <- apply(Yadj, 2, function(col) cor(s, col))
  expect_true(mean(abs(cors), na.rm = TRUE) < 0.05)
})

test_that("pc_remove: study mode skips small studies and stitches back in order", {
  set.seed(6)
  G <- 10
  # Two studies: A has 12 samples; B has 4 samples (below default min_samples=15, we'll set 6)
  A <- 12; B <- 4
  S <- A + B
  M <- matrix(rnorm(G * S), nrow = G, ncol = S)
  rownames(M) <- paste0("g", seq_len(G))
  colnames(M) <- c(paste0("A", 1:A), paste0("B", 1:B))

  study_map <- data.frame(
    Experiment = colnames(M),
    Study      = c(rep("A", A), rep("B", B)),
    stringsAsFactors = FALSE
  )

  # Set min_samples so A is corrected and B is skipped
  out <- pc_remove(M, study_map = study_map, method = "study",
                   n_pcs = 1L, min_samples = 6L)

  X <- out$te_corrected

  # Columns are stitched back in original order
  expect_identical(colnames(X), colnames(M))

  # For the skipped study B, the block must be unchanged
  expect_equal(X[, paste0("B", 1:B), drop = FALSE],
               M[, paste0("B", 1:B), drop = FALSE])

  # meta table present with A having n_pcs >= 0 and B == 0
  expect_true(is.data.frame(out$n_pcs))
  expect_setequal(out$n_pcs$Study, c("A", "B"))
  expect_identical(out$n_pcs$n_pcs[out$n_pcs$Study == "B"], 0L)
})

test_that("pc_remove: dimensions and names always preserved", {
  set.seed(7)
  G <- 9; S <- 11
  M <- matrix(rnorm(G * S), nrow = G, ncol = S)
  rownames(M) <- paste0("gene_", seq_len(G))
  colnames(M) <- paste0("sample_", seq_len(S))

  out <- pc_remove(M, method = "global", n_pcs = 2L)
  X   <- out$te_corrected

  expect_identical(dim(X), dim(M))
  expect_identical(rownames(X), rownames(M))
  expect_identical(colnames(X), colnames(M))
})

# tests/testthat/test-pc_remove_edges.R

# -- Helpers ---------------------------------------------------------------

proj_residuals <- function(Y, X) {
  # Residuals: Y - X (X'X)^(-1) X' Y
  beta <- solve(crossprod(X), crossprod(X, Y))
  Y - X %*% beta
}

# -- Input validation edge cases ------------------------------------------

test_that("errors: empty, single-dim, non-finite, bad map", {
  # empty dims
  M0 <- matrix(numeric(0), nrow = 0, ncol = 0,
               dimnames = list(character(0), character(0)))
  expect_error(pc_remove(M0, method = "global"), "rownames|colnames", ignore.case = TRUE)

  # NA/Inf
  M3 <- matrix(rnorm(6), nrow = 2)
  rownames(M3) <- c("g1","g2"); colnames(M3) <- c("s1","s2","s3")
  M3[1,1] <- NA_real_
  expect_error(pc_remove(M3, method = "global"), "NA/Inf", ignore.case = TRUE)

  # bad study_map
  M4 <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(M4) <- paste0("g", 1:3); colnames(M4) <- paste0("s", 1:4)
  bad_map <- data.frame(Experiment = paste0("s", 1:3), Study = "A")
  expect_error(pc_remove(M4, method = "study", study_map = bad_map),
               "Missing.*study_map\\$Experiment", ignore.case = TRUE)
})

# -- Boundary k behavior ---------------------------------------------------

test_that("k < 0, NA, or huge are handled by clamping", {
  set.seed(10)
  G <- 7; S <- 5
  M <- matrix(rnorm(G*S), nrow = G,
              dimnames = list(paste0("g",1:G), paste0("s",1:S)))

  out1 <- pc_remove(M, method = "global", n_pcs = -5L)
  expect_identical(out1$n_pcs, 0L)

  out2 <- pc_remove(M, method = "global", n_pcs = as.integer(NA))
  expect_identical(out2$n_pcs, 0L)

  out3 <- pc_remove(M, method = "global", n_pcs = 9999L)
  expect_identical(out3$n_pcs, as.integer(min(S, G) - 1L))
})

# -- Orthogonality property: residuals ⟂ removed subspace ------------------

test_that("residuals are orthogonal to the removed U_k", {
  set.seed(11)
  G <- 20; S <- 10
  M <- matrix(rnorm(G*S), nrow = G,
              dimnames = list(paste0("g",1:G), paste0("s",1:S)))

  # Correct with k=2 and reconstruct the U used internally to verify orthogonality.
  # We replicate the steps here (not touching pc_remove internals).
  # NOTE: zero-variance isn’t present in random data.
  Y <- t(M)                    # samples x genes
  sv <- svd(Y, nu = 2, nv = 0)
  U  <- sv$u[, 1:2, drop = FALSE]
  X  <- cbind(1, U)
  Y_adj_expected <- proj_residuals(Y, X)

  # Call our function
  out <- pc_remove(M, method = "global", n_pcs = 2L)
  Y_adj <- t(out$te_corrected)

  # Orthogonality: U' * Y_adj ≈ 0  (each removed PC has ~0 correlation with residuals)
  Z <- t(U) %*% Y_adj     # 2 x genes
  expect_true(max(abs(Z)) < 1e-8)
  # And the residuals match the explicit projector (up to numerical tolerance)
  expect_equal(Y_adj, Y_adj_expected, tolerance = 1e-8)
})

# -- Idempotence: applying k=0 twice does nothing further ------------------

test_that("idempotence for k=0", {
  set.seed(12)
  G <- 12; S <- 7
  M <- matrix(rnorm(G*S), nrow = G,
              dimnames = list(paste0("g",1:G), paste0("s",1:S)))

  out1 <- pc_remove(M, method = "global", n_pcs = 0L)$te_corrected
  out2 <- pc_remove(out1, method = "global", n_pcs = 0L)$te_corrected
  expect_equal(out1, out2, tolerance = 1e-12)  # exactly equal
})

# -- Zero-variance handling: all-constant block ----------------------------

test_that("all-constant matrix is returned unchanged and k=0", {
  G <- 5; S <- 4
  M <- matrix(7, nrow = G, ncol = S,
              dimnames = list(paste0("g",1:G), paste0("s",1:S)))

  out <- pc_remove(M, method = "global", n_pcs = 1L)
  expect_identical(out$te_corrected, M)
  expect_identical(out$n_pcs, 0L)  # nothing to remove
  expect_setequal(out$dropped_genes, rownames(M))  # all were dropped for SVD
})

# -- Order preservation under shuffles ------------------------------------

test_that("row/col order is preserved", {
  set.seed(13)
  G <- 9; S <- 8
  M <- matrix(rnorm(G*S), nrow = G,
              dimnames = list(paste0("g",1:G), paste0("s",1:S)))

  # Shuffle rows/cols
  rperm <- sample(1:G); cperm <- sample(1:S)
  M2 <- M[rperm, cperm, drop = FALSE]

  out <- pc_remove(M2, method = "global", n_pcs = 1L)$te_corrected

  # Output must match the shuffled ordering, not revert to original
  expect_identical(rownames(out), rownames(M2))
  expect_identical(colnames(out), colnames(M2))
})

# -- Study mode: overlapping/misaligned labels and boundary of min_samples --

test_that("study mode respects min_samples boundary and alignment", {
  set.seed(14)
  G <- 6
  A <- 6  # exactly equals min_samples
  B <- 5  # below min_samples
  S <- A + B

  M <- matrix(rnorm(G*S), nrow = G,
              dimnames = list(paste0("g",1:G),
                              c(paste0("A",1:A), paste0("B",1:B))))

  study_map <- data.frame(
    Experiment = colnames(M),
    Study      = c(rep("A", A), rep("B", B)),
    stringsAsFactors = FALSE
  )

  out <- pc_remove(M, method = "study", study_map = study_map,
                   n_pcs = 1L, min_samples = 6L)

  # A corrected (>= boundary), B unchanged (< boundary)
  expect_true(is.data.frame(out$n_pcs))
  expect_identical(sort(out$n_pcs$Study), c("A","B"))
  expect_identical(out$n_pcs$n_pcs[out$n_pcs$Study=="B"], 0L)
  expect_false(identical(out$te_corrected[, paste0("A",1:A)], M[, paste0("A",1:A)]))
  expect_identical(out$te_corrected[, paste0("B",1:B)], M[, paste0("B",1:B)])
})

# -- k equals max possible: residual rank sanity ---------------------------

test_that("k = max yields only intercept residuals (demeaned samples)", {
  set.seed(15)
  G <- 10; S <- 4
  M <- matrix(rnorm(G*S), nrow = G,
              dimnames = list(paste0("g",1:G), paste0("s",1:S)))

  k_max <- as.integer(min(S, G) - 1L)  # = 3
  out <- pc_remove(M, method = "global", n_pcs = k_max)
  X   <- out$te_corrected
  # After removing all non-intercept subspace, each sample is demeaned:
  expect_equal(as.numeric(colMeans(X)), rep(0, S), tolerance = 1e-10)
})

# -- Auto mode (optional): skip if sva missing -----------------------------

test_that("n_pcs = 'auto' behaves (if sva installed)", {
  testthat::skip_if_not_installed("sva")
  set.seed(16)
  G <- 12; S <- 10
  M <- matrix(rnorm(G*S), nrow = G,
              dimnames = list(paste0("g",1:G), paste0("s",1:S)))
  out <- pc_remove(M, method = "global", n_pcs = "auto")
  expect_true(is.integer(out$n_pcs))
  expect_gte(out$n_pcs, 0L)
  expect_lte(out$n_pcs, as.integer(min(S, G) - 1L))
})
