#' Select subset of genes for downstream processing
#'
#' Assumptions
#' 1) Inputs are genes by samples
#' 2) Both RNA and Ribo are calculated for same set of samples
#'
#' Existing features
#' a) If `dummy = TRUE` the counts of removed genes are summed into a new row
#' b) CPM is computed column-wise in each matrix and evaluated over the combined columns
#'
#' Optional features to Implement later
#' a) Support only ribo or RNA
#' b) Support removing genes without polyA genes
#'
#' @param RIBO numeric matrix or data.frame with genes in rows and samples in columns
#' @param RNA numeric matrix or data.frame with genes in rows and samples in columns
#' @param dummy logical. If TRUE add a row that contains the column sums of removed genes
#' @param cpm numeric CPM threshold. Default 1
#' @param fraction  numeric between 0,1. A gene is kept only if at least this fraction
#' of the combined columns has CPM >= `cpm`. Default 0.8
#'
#' @return list with elements:
#'   - ribo: filtered RIBO matrix (genes by samples, with optional dummy row)
#'   - rna:  filtered RNA matrix (genes by samples, with optional dummy row)
#'   - kept_genes: character vector of kept gene ids
#'   - removed_genes: character vector of removed gene ids
#' @examples
#' set.seed(3)
#' G <- 100; S <- 8
#' RIBO <- matrix(rpois(G*S, 20), nrow = G,
#'                dimnames = list(paste0("g",1:G), paste0("s",1:S)))
#' RNA  <- matrix(rpois(G*S, 25), nrow = G,
#'                dimnames = list(paste0("g",1:G), paste0("s",1:S)))
#' # make some low CPM genes
#' RIBO[1:5,] <- 0
#' RNA [1:5,] <- 0
#' sel <- select_genes(RIBO, RNA, cpm = 1, fraction = 0.8, dummy = TRUE)
#' dim(sel$ribo); dim(sel$rna)
#' @importFrom stats setNames
#' @export

select_genes <- function(RIBO, RNA,
                                 dummy = TRUE,
                                 cpm = 1,
                                 fraction = 0.8) {
  RIBO <- as.matrix(RIBO)
  RNA  <- as.matrix(RNA)

  # intersect genes and samples to guarantee alignment
  common_genes  <- intersect(rownames(RIBO), rownames(RNA))
  if (length(common_genes) == 0L) stop("No shared genes between RIBO and RNA")
  common_samples <- intersect(colnames(RIBO), colnames(RNA))
  if (length(common_samples) == 0L) stop("No shared samples between RIBO and RNA")

  RIBO <- RIBO[common_genes, common_samples, drop = FALSE]
  RNA  <- RNA [common_genes, common_samples, drop = FALSE]

  cpm_fun <- function(M) {
    lib <- colSums(M)
    sweep(M, 2L, lib, "/") * 1e6
  }

  # CPM per modality
  ribo_cpm <- cpm_fun(RIBO)
  rna_cpm  <- cpm_fun(RNA)

  # evaluate expression over the combined columns
  combined_cpm <- cbind(ribo_cpm, rna_cpm)  # genes by (samples_RIBO + samples_RNA)
  # keep gene if at least 'fraction' of combined columns are >= cpm
  prop_ge <- rowMeans(combined_cpm >= cpm)
  keep <- prop_ge >= fraction

  kept_genes    <- rownames(RIBO)[keep]
  removed_genes <- rownames(RIBO)[!keep]

  RIBO_keep <- RIBO[keep, , drop = FALSE]
  RNA_keep  <- RNA [keep, , drop = FALSE]

  if (isTRUE(dummy) && length(removed_genes) > 0L) {
    # collapse removed counts into a single row for each matrix
    dummy_ribo <- colSums(RIBO[!keep, , drop = FALSE])
    dummy_rna  <- colSums(RNA [!keep, , drop = FALSE])
    RIBO_keep  <- rbind(RIBO_keep, DUMMY_REMOVED = dummy_ribo)
    RNA_keep   <- rbind(RNA_keep,  DUMMY_REMOVED = dummy_rna)
  }

  list(
    ribo = RIBO_keep,
    rna  = RNA_keep,
    kept_genes = kept_genes,
    removed_genes = removed_genes
  )

}
