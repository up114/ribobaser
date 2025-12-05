#' Principal-component removal for TE matrices (global or study-level)
#'
#' Remove leading principal components (PCs) from a translation efficiency (TE)
#' matrix to control unwanted sample-level structure. Works either globally
#' (one SVD across all samples) or per study (one SVD per study block).
#'
#' Assumptions
#'  1) Inputs are genes by samples
#'  2) Inputs are already CLR-normalized (no normalization occurs within this function)
#'    - either run with TE matrix or with Ribo-Seq/RNA-Seq matrices after CLR normalization
#'
#' @param mat numeric matrix, genes x samples. Row names = genes; col names = samples.
#' @param study_map optional data.frame with columns {Experiment} and {Study};
#'   required when {method = "study"}. {Experiment} must match {colnames(mat)}.
#' @param method character, {"study"} or {"global"}.
#' @param n_pcs {"auto"} (default), which uses sva::num.sv(..., method = "be"), or non-negative integer.
#' @param min_samples integer; for {method="study"}, skip correction when a study
#'   has fewer than this many samples. Default: 15.
#'
#' @return list with elements:
#'    - mat_pc: matrix (genes x samples), with the same dimensions as the input matrix
#'    - n_pcs: if global, an integer; if study, a data.frame with {Study} and {n_pcs}.
#'    - dropped_genes: character vector of genes dropped from SVD due to near-zero variance
#'         (they are re-inserted unchanged in the output).
#'
#' @examples
#' # Global correction:
#' # out <- pc_remove(mat, method = "global")
#'
#' @importFrom stats lm residuals
#' @export
pc_remove <- function(te_mat,
                      study_map = NULL,
                      method = c("study", "global"),
                      n_pcs = "auto",
                      min_samples = 15) {

  method <- match.arg(method)

  # ----- checks -----
  stopifnot(is.matrix(te_mat), is.numeric(te_mat))
  if (is.null(rownames(te_mat)) || is.null(colnames(te_mat))) {
    stop("te_mat must have rownames (genes) and colnames (samples).")
  }
  if (method == "study") {
    if (is.null(study_map) || !all(c("Experiment","Study") %in% colnames(study_map))) {
      stop("Provide study_map with columns: Experiment, Study.")
    }
    study_map$Experiment <- as.character(study_map$Experiment)
    study_map$Study      <- as.character(study_map$Study)
    if (!all(colnames(te_mat) %in% study_map$Experiment)) {
      miss <- setdiff(colnames(te_mat), study_map$Experiment)
      stop("Missing in study_map$Experiment: ", paste(miss, collapse = ", "))
    }
  }
  if (any(!is.finite(te_mat))) stop("te_mat contains NA/Inf; handle before pc_remove().")

  # -----  helper functions -----

  # Drop ~zero-variance genes (for SVD stability). They’re re-inserted unchanged later.
  drop_zero_var <- function(M, tol = 1e-6) {
    v <- matrixStats::rowVars(M)

    keep <- !(v <= tol | is.na(v))
    list(M = M[keep, , drop = FALSE], keep = keep)
  }


  # Decide how many PCs to remove:
  choose_k <- function(Y) {
    # Y is samples x genes - transpose
    if (identical(n_pcs, "auto")) {
      if (!requireNamespace("sva", quietly = TRUE)) {
        stop("Install 'sva' for n_pcs='auto'.")
      }
      mod <- matrix(1, nrow = nrow(Y), ncol = 1)
      k   <- as.integer(sva::num.sv(t(Y), mod, method = "be"))
    } else {
      k <- as.integer(n_pcs)
    }
    if (is.na(k) || k < 0) k <- 0L
    # clamp: keep k within [0, min(dim(Y)) - 1]
    max(0L, min(k, min(dim(Y)) - 1L))
  }

  # Correct one block (either the whole matrix, or one study subset)
  correct_block <- function(block, label = NULL) {
    cat("== Block:", ifelse(is.null(label), "GLOBAL", label),
        "| dims:", nrow(block), "genes x", ncol(block), "samples\n")

    # 1) drop ~zero-variance genes
    dz <- drop_zero_var(block)
    M  <- dz$M
    kept_idx <- dz$keep
    cat("   dropped zero-variance genes:", sum(!kept_idx), "\n")

    if (nrow(M) == 0) {
      cat("   all genes ~zero-variance; skipping\n")
      return(list(corrected = block, k = 0L, dropped = rownames(block)[!kept_idx]))
    }

    # 2) transpose: samples x genes
    Y <- t(M)

    # 3) pick k
    k <- choose_k(Y)
    cat("   PCs to remove:", k, "\n")

    # 4) regress out top k PCs (equivalent to lm(Y ~ U) with intercept)
    if (k > 0) {
      sv  <- svd(Y, nu = k, nv = 0)               # only need U (left singular vectors)
      U   <- sv$u[, seq_len(k), drop = FALSE]
      X   <- cbind(1, U)                          # design = [Intercept, U_k]
      beta <- solve(crossprod(X), crossprod(X, Y))
      Y_adj <- Y - X %*% beta                     # residuals
    } else {
      # no PCs: just demean each sample (row)
      Y_adj <- Y - rowMeans(Y)
    }

    # 5) back to genes x samples; re-insert dropped genes unchanged
    adj <- t(Y_adj)
    out <- block
    out[kept_idx, ] <- adj

    list(corrected = out, k = k, dropped = rownames(block)[!kept_idx])
  }

  # ----- main flow -----
  original_cols <- colnames(te_mat)
  dropped_all   <- character(0)

  if (method == "global") {
    g <- correct_block(te_mat, "GLOBAL")
    te_corr <- g$corrected
    n_out   <- g$k
    dropped_all <- g$dropped

  } else {
    # study-wise correction
    map <- setNames(study_map$Study, study_map$Experiment)
    studies <- unique(map[original_cols])                 # preserve input order
    blocks  <- vector("list", length(studies))
    names(blocks) <- studies
    meta <- data.frame(Study = studies, n_pcs = NA_integer_, stringsAsFactors = FALSE)

    for (s in studies) {
      cols <- names(which(map == s))
      sub  <- te_mat[, cols, drop = FALSE]
      if (ncol(sub) < min_samples) {
        cat("== Block:", s, "| too few samples (", ncol(sub),
            " < ", min_samples, ") -> skipping\n", sep = "")
        blocks[[s]] <- sub
        meta$n_pcs[meta$Study == s] <- 0L
      } else {
        res <- correct_block(sub, s)
        blocks[[s]] <- res$corrected
        meta$n_pcs[meta$Study == s] <- res$k
        dropped_all <- c(dropped_all, res$dropped)
      }
    }

    te_corr <- do.call(cbind, blocks)[, original_cols, drop = FALSE]
    n_out   <- meta
  }

  # match input ordering
  te_corr <- te_corr[rownames(te_mat), original_cols, drop = FALSE]

  list(
    te_corrected  = te_corr,
    n_pcs         = n_out,
    dropped_genes = unique(dropped_all)
  )
}
