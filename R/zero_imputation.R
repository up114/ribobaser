#' Replace zero counts with pseudocounts
#'
#' Many compositional workflows require strictly positive counts before
#' applying log-ratio transforms. This utility performs a simple zero-imputation
#' step, replacing all zeros in paired RIBO and RNA matrices with ones.
#'
#' Features to implement
#' a) Multiplicative imputation as in Shilpa's work
#'
#' @param RIBO numeric matrix or data.frame of counts with genes in rows and
#'   samples in columns.
#' @param RNA numeric matrix or data.frame of counts with genes in rows and
#'   samples in columns.
#'
#' @return list with elements
#'   - ribo: `RIBO` converted to a matrix with zeros replaced by ones.
#'   - rna:  `RNA` converted to a matrix with zeros replaced by ones.
#'
#' @examples
#' RIBO <- matrix(c(0, 5, 10,
#'                  3, 0,  7),
#'                nrow = 2,
#'                dimnames = list(c("g1", "g2"), paste0("s", 1:3)))
#' RNA  <- matrix(c(2, 0, 4,
#'                  0, 6, 8),
#'                nrow = 2,
#'                dimnames = list(c("g1", "g2"), paste0("s", 1:3)))
#'
#' zero_imputation(RIBO, RNA)
#'
#' @export
zero_imputation <- function(RIBO, RNA) {
  ribo_mat <- as.matrix(RIBO)
  rna_mat  <- as.matrix(RNA)

  ribo_mat[ribo_mat == 0] <- 1
  rna_mat [rna_mat  == 0] <- 1

  list(
    ribo = ribo_mat,
    rna  = rna_mat
  )
}
