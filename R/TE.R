#' Translation efficiency via compositional regression or CLR difference
#'
#' Compute TE from paired ribosome profiling and
#' RNA-seq counts under compositional data analysis principles.
#'
#' This function supports two  approaches controlled by `method`:
#' * `'regression'` (default) reproduces the Liu et al 2025 workflow
#' that converts CLR values to ILR coordinates and
#' regresses ILR(RIBO) on ILR(RNA) per sample, mapping residuals back to CLR space.
#' * `'logratio'` computes the simpler CLR difference `clr(RIBO) - clr(RNA)`
#'
#' Assumptions
#' 1) Inputs are samples by genes.
#' 2) Zeros were handled before calling this function.
#' 3) No internal scaling or centering is performed.
#'
#' Things to Check
#' 1) Consider switching out propr call with the simple function. Make sure the equivalence/zero handling
#'
#' @param RIBO numeric matrix or data.frame with samples in rows and genes in columns
#' @param RNA numeric matrix or data.frame with samples in rows and genes in columns
#' @param method character string selecting the TE estimator. `'regression'`
#'   runs the compositional regression workflow, `'logratio'` returns
#'   `clr(RIBO) - clr(RNA)`.
#' @param pc_removal logical for whether principal component removal should be applied
#'        to the RIBO and RNA matrices after conversion to CLR values
#' @param n_pcs {"auto"} (default), which uses sva::num.sv(..., method = "be"), or non-negative integer.
#' @param parallel logical run with foreach dopar (only used for
#'   `method = "regression"`).
#' @param n_cores integer cores for parallel (only used for
#'   `method = "regression"`).
#'
#' @return numeric matrix genes by samples
#' @examples
#' RIBO <- matrix(c(5, 7, 9,
#'                  6, 8, 10),
#'                nrow = 2,
#'                byrow = TRUE,
#'                dimnames = list(paste0("g", 1:2), paste0("s", 1:3)))
#' RNA <- matrix(c(4, 6, 8,
#'                 5, 7, 9),
#'               nrow = 2,
#'               byrow = TRUE,
#'               dimnames = list(paste0("g", 1:2), paste0("s", 1:3)))
#' te(RIBO, RNA, method = "logratio")
#'
#'
#' @importFrom stats lm residuals
#' @importFrom foreach %do% %dopar%
#' @export
te <- function(RIBO, RNA,
                  method = c("regression", "logratio"),
                  pc_removal = FALSE,
                  n_pcs = "auto",
                  parallel = FALSE,
                  n_cores = max(1L, parallel::detectCores() - 1L)) {

  method <- match.arg(method)

  RIBO <- as.matrix(RIBO)
  RNA  <- as.matrix(RNA)

  # intersect samples and genes
  if (!identical(rownames(RIBO), rownames(RNA))) {
    common_samples <- intersect(rownames(RIBO), rownames(RNA))
    if (length(common_samples) == 0L) stop("No shared samples")
    RIBO <- RIBO[common_samples, , drop = FALSE]
    RNA  <- RNA [common_samples, , drop = FALSE]
  }
  if (!identical(colnames(RIBO), colnames(RNA))) {
    common_genes <- intersect(colnames(RIBO), colnames(RNA))
    if (length(common_genes) == 0L) stop("No shared genes")
    RIBO <- RIBO[, common_genes, drop = FALSE]
    RNA  <- RNA [, common_genes, drop = FALSE]
  }

  if (method == "logratio") {
    if (any(RIBO <= 0) || any(RNA <= 0)) {
      stop("Method 'logratio' requires strictly positive counts; handle zeros before calling")
    }
    clr_transform <- function(mat) {
      log_mat <- log(mat)
      sweep(log_mat, 1, rowMeans(log_mat), "-")
    }
    te_clr <- t(clr_transform(RIBO) - clr_transform(RNA))
    colnames(te_clr) <- rownames(RIBO)
    rownames(te_clr) <- colnames(RIBO)
    return(as.matrix(te_clr))
  }

  # CLR via propr
  pr_RIBO <- propr::propr(RIBO, metric = "rho", ivar = "clr", alpha = NA, p = 100)
  pr_RNA  <- propr::propr(RNA,  metric = "rho", ivar = "clr", alpha = NA, p = 100)

  if(pc_removal) {
    ribo_clr_pc <- pc_remove(t(pr_RIBO@logratio),
                             study_map = NULL,
                             method = "global",
                             n_pcs = n_pcs,
                             min_samples = 15)

    rna_clr_pc <- pc_remove(t(pr_RNA@logratio),
                            study_map = NULL,
                            method = "global",
                            n_pcs = n_pcs,
                            min_samples = 15)
    ribo_clr_pc <- t(ribo_clr_pc$te_corrected)
    rna_clr_pc <- t(rna_clr_pc$te_corrected)

    # force row means back to 0
    pr_RIBO@logratio <- sweep(ribo_clr_pc, 1, rowMeans(ribo_clr_pc), "-")
    pr_RNA@logratio <- sweep(rna_clr_pc, 1, rowMeans(rna_clr_pc), "-")
  }

  # CLR -> ILR
  RIBO_ilr <- compositions::clr2ilr(pr_RIBO@logratio)
  RNA_ilr  <- compositions::clr2ilr(pr_RNA@logratio)

  # transposed to ILRcoords by samples like your script
  RIBO_ilr <- as.data.frame(t(RIBO_ilr))
  RNA_ilr  <- as.data.frame(t(RNA_ilr))

  # parallel setup
  `%op%` <- if (isTRUE(parallel)) `%dopar%` else `%do%`
  if (isTRUE(parallel)) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
    }, add = TRUE)
  }

  # per sample regression across ILR coordinates
  out <- foreach::foreach(i = 1:ncol(RIBO_ilr),
                          .combine = "cbind",
                          .packages = c("compositions", "stats"),
                          .inorder = TRUE) %op% {
                            m <- stats::lm(RIBO_ilr[, i] ~ RNA_ilr[, i])
                            # residuals live in ILR space over coordinates
                            r_ilr <- stats::residuals(m)
                            # back to CLR per sample
                            r_clr <- as.numeric(compositions::ilr2clr(r_ilr))
                            data.frame(r_clr)
                          }

  # names to match
  colnames(out) <- rownames(RIBO)  # samples as columns
  rownames(out) <- colnames(RIBO)  # genes as rows

  as.matrix(out)
}
