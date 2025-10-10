#' Translation efficiency covariation
#'
#' Compute translation efficiency covariation (TEC) for all gene pairs from a matrix of
#' translation efficiency (TE) values. Input should be genes in rows and samples in columns.
#' The current implementation relies on `propr:::lr2rho()` to obtain proportionality scores.
#'
#' The computation scales quadratically with the number of genes because every pair is
#' evaluated. For large gene sets and parallelization strategies.
#'
#' Assumptions
#' 1) `TE` is obtained via `te(RIBO, RNA, method = "regression")`, i.e. values are residual
#'    CLR coordinates in gene space.
#' 2) Rows correspond to genes, columns to samples that already share a common coordinate
#'    system.
#' 3) No missing values remain after earlier preprocessing steps.
#'
#' @param TE numeric matrix or data.frame with genes in rows and samples in columns.
#' @param method character string selecting the TEC estimator. Currently only `'rho'`
#'   is available and it uses the `propr` package.
#'
#' @return Symmetric numeric matrix genes by genes with proportionality scores.
#' @examples
#'   data("ribo_raw_human_cap_995", package = "ribobaser")
#'   data("rnaseq_raw_human_cap_995", package = "ribobaser")
#'
#'   keep_genes <- rownames(ribo_raw_human_cap_995)[4:7]
#'   keep_samples <- colnames(ribo_raw_human_cap_995)[1:5]
#'
#'   ribo <- t(as.matrix(ribo_raw_human_cap_995[keep_genes, keep_samples]))
#'   rna  <- t(as.matrix(rnaseq_raw_human_cap_995[keep_genes, keep_samples]))
#'
#'   te_mat <- te(ribo, rna, method = "regression")
#'   tec(te_mat)
#' 
#'
#' @export
tec <- function(TE,
                method = c("rho")) {

  method <- match.arg(method)

  if (!requireNamespace("propr", quietly = TRUE)) {
    stop("Package 'propr' is required for method 'rho'. Install it with install.packages('propr').",
         call. = FALSE)
  }

  TE <- as.matrix(TE)

  if (anyNA(TE)) {
    stop("TEC requires TE values without missing data; impute or remove NAs before calling tec().")
  }
  if (!all(is.finite(TE))) {
    stop("TEC detected non-finite values in TE; ensure upstream TE computation removed infinities.")
  }

  gene_count <- nrow(TE)
  if (gene_count < 2L) {
    stop("TEC requires at least two genes (rows) in the TE matrix.")
  }

  gene_names <- rownames(TE)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene", seq_len(gene_count))
  }

  lr2rho <- getFromNamespace("lr2rho", "propr")
  rho_vals <- lr2rho(t(TE))

  if (is.matrix(rho_vals) && nrow(rho_vals) == gene_count && ncol(rho_vals) == gene_count) {
    rho_mat <- rho_vals
  } else {
    pair_count <- gene_count * (gene_count - 1L) / 2L
    if (length(rho_vals) != pair_count) {
      stop("Unexpected output from propr:::lr2rho; expected vector of length choose(n_genes, 2).")
    }
    rho_mat <- matrix(1, nrow = gene_count, ncol = gene_count)
    rho_mat[upper.tri(rho_mat)] <- rho_vals
    rho_mat[lower.tri(rho_mat)] <- t(rho_mat)[lower.tri(rho_mat)]
  }

  dimnames(rho_mat) <- list(gene_names, gene_names)
  rho_mat
}
