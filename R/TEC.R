#' Translation efficiency covariation calculation
#'
#' Compute TEC given a matrix of TE values by samples.
#' Samples may have been aggregated by cell type or any other criteria.
#' Note that this computation can be computationally intensive.
#' Assumptions
#' 1) Input is genes by samples
#' 2) TE calculation for each sample was done and samples were aggregated
#'
#' Steps
#' a) Yet to be implemented
#'
#' Features to Implement
#' a) Add alternative measures of TEC.
#' b) Add parallelization support
#' c) Better formatting of the return value. Consider just the upper diag
#' @param TE numeric matrix or data.frame with genes in rows and samples in columns
#'
#' @return numeric matrix genes by genes. Each value is TEC between the gene pairs.
#' @examples
#' # Example with simple TE matrix
#'
#' @export
