#' Translation efficiency covariation
#'
#' Compute translation efficiency covariation (TEC) for all gene pairs from a matrix of
#' translation efficiency (TE) values. Input should be genes in rows and samples in columns.
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
#' @param method character string selecting the TEC estimator.
#'   - `'rho'` uses the `propr` package and returns a TEC matrix;
#'     `propr:::lr2rho()` to obtain proportionality scores
#'   - `'wgcna'` uses the `WGCNA` package and returns a TOM matrix
#'   - `'glasso'` uses the `huge` package and returns a precision matrix
#'   - `'genie3'` uses the `GENIE3` package and returns a weight matrix
#' @param n_cores integer cores for parallel (only used for
#'   `method = "genie3"`).
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
#'   tec_mat = tec(te_mat, method = "rho")
#'
#' @export
tec <- function(TE,
                method = c("rho", "glasso", "wgcna", "genie3"),
                n_cores = max(1L, parallel::detectCores() - 1L)) {

  method <- match.arg(method)

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

  if(method == "rho") {

    if (!requireNamespace("propr", quietly = TRUE)) {
      stop("Package 'propr' is required for method 'rho'. Install it with install.packages('propr').",
           call. = FALSE)
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
    return(rho_mat)
  }

  else if(method == "glasso") {
    if (!requireNamespace("huge", quietly = TRUE)) {
      stop("Package 'huge' is required for method 'glasso'. Install it with install.packages('huge').",
           call. = FALSE)
    }

    TE <- t(TE)
    keep <- apply(TE, 2, sd, na.rm = TRUE) > 1e-6
    if (!any(keep)) {
      stop("All genes have near-zero variance; cannot estimate a glasso TEC network.")
    }

    TE <- TE[, keep, drop = FALSE]

    # center and scale each gene for numerical stability
    TE <- scale(TE)

    lambda_grid = seq(0.1, 1.00, by = 0.02)


    fit = huge::huge(TE,
                     lambda = lambda_grid,
                     method = "glasso",
                     scr = TRUE,
                     verbose = FALSE
    )

    opt_graph = huge::huge.select(
      fit,
      criterion = "ric",
      verbose = FALSE
    )

    lambda <- opt_graph$opt.lambda
    cat(sprintf("Chosen lambda = %.2f\n", lambda))

    Theta = opt_graph$opt.icov # precision matrix

    # Save results
    Theta <- (Theta + t(Theta)) / 2 # symmetrize

    gene_names <- colnames(TE)
    if (is.null(gene_names)) {
      gene_names <- paste0("gene", seq_len(ncol(TE)))
    }
    dimnames(Theta) <- list(gene_names, gene_names)

    return(Theta)
  }

  else if(method == "wgcna") {
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      stop("Package 'WGCNA' is required for method 'wgcna'. Install it with install.packages('WGCNA').",
           call. = FALSE)
    }
    TE <- t(TE)

    # Basic QC for NAs/zeros
    # removes genes/samples with too many missing values or zero variance
    gsg <- WGCNA::goodSamplesGenes(TE, verbose = 0)
    if (!gsg$allOK) {
      TE <- TE[gsg$goodSamples, gsg$goodGenes]
    }

    if(ncol(TE) < 1) {
      stop("All genes dropped in WGCNA QC.")
    }

    # Soft-threshold (power) scan
    powers <- c(1:10, seq(12, 20, 2))  # candidate powers for adjacency function

    sft <- WGCNA::pickSoftThreshold(
      TE,
      networkType = "signed",                     # signed network
      corFnc = "cor",                             # correlation function used to compute similarity
      corOptions = list(use = "pairwise.complete.obs"),  # handles missing data by pairwise deletion
      powerVector = powers,                       # the vector of powers to test
      verbose = 0
    )

    # pick power (find minimum power where R^2 >= 0.9)
    ok <- which(sft$fitIndices$SFT.R.sq >= 0.9)
    if (length(ok) == 0) {
      warning("No power reaches R^2 >= 0.9; picking power with max R^2.")
      picked_power <- sft$fitIndices$Power[ which.max(sft$fitIndices$SFT.R.sq) ]
    } else {
      picked_power <- sft$fitIndices$Power[ ok[1] ]
    }

    cat(sprintf("Chosen power = %.2f\n", picked_power))

    # generate correlation matrix - TOM
    TOM <- WGCNA::TOMsimilarityFromExpr(
      TE,
      networkType = "signed",
      TOMType = "signed",
      power = picked_power,
      corType = "pearson"
    )

    # Format + return genes×genes matrix
    gene_names <- colnames(TE)
    if (is.null(gene_names)) {
      gene_names <- paste0("gene", seq_len(ncol(TE)))
    }

    dimnames(TOM) <- list(gene_names, gene_names)
    TOM <- (TOM + t(TOM)) / 2
    return(TOM)
  }
  else if(method == "genie3") {
    if (!requireNamespace("GENIE3", quietly = TRUE)) {
      stop("Package 'GENIE3' is required for method 'genie3'. Install it with BiocManager::install('GENIE3').",
           call. = FALSE)
    }

    set.seed(42)
    W <- GENIE3::GENIE3(as.matrix(TE), nCores = n_cores, nTrees = 100)

    # Check same gene set (ignoring order)
    if (!setequal(rownames(W), gene_names)) {
      stop("GENIE3 output rows do not match the TE gene set.")
    }
    if (!setequal(colnames(W), gene_names)) {
      stop("GENIE3 output columns do not match the TE gene set.")
    }

    #  reordering to original TE order
    W <- W[gene_names, gene_names, drop = FALSE]

    # Make it symmetric by averaging both directions
    W_sym <- (W + t(W)) / 2
    dimnames(W_sym) <- list(gene_names, gene_names)

    return(W_sym)
  }
}
