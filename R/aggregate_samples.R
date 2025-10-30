#' Aggregate Samples
#'
#' Samples can be aggregated by cell line or study using the group mean.
#' Optional features to Implement later
#' a) Support only ribo or RNA
#' b) Provide option to use either Ribobase_QC_non_dedup_data or Ribobase_QC_dedup_data
#' c) Provide support for other mathematical functions to aggregate the numerical values
#' d) Need to explore Battle_GR paper for additional ideas.
#' e) Potentially add SVA/PCA based confounder correction. This might be in a different script
#'
#' Given TE matrix with genes in rows and samples in
#' columns, collapse samples according to metadata supplied in
#' `Ribobase_QC_dedup_data`. The function currently supports aggregating by
#' `Study` (GSE accession) or by `Cell line`, computing the mean TE for each gene
#' within the chosen grouping.
#'
#' @param TE numeric matrix or data.frame with genes in rows and samples in
#'   columns. Column names must match the `Experiment` identifiers present in
#'   `Ribobase_QC_dedup_data`.
#' @param group_by character string specifying the metadata field to use when
#'   combining samples. Either `'study'` or `'cell_line'`. Defaults to `'study'`.
#' @param fun aggregation function identifier. Currently only `'mean'` is
#'   supported.
#'
#' @return A numeric matrix with genes in rows and aggregated sample groups in
#'   columns. Column names correspond to the selected metadata grouping.
#'
#' @examples
#' data("Ribobase_QC_dedup_data", package = "ribobaser")
#' example_samples <- Ribobase_QC_dedup_data$Experiment[1:4]
#' TE <- matrix(
#'   seq_along(example_samples) + rep(0:1, each = length(example_samples)),
#'   nrow = 2,
#'   dimnames = list(paste0("g", 1:2), example_samples)
#' )
#' aggregate_samples(TE, group_by = "study")
#' aggregate_samples(TE, group_by = "cell_line")
#' @importFrom stats setNames
#' @importFrom utils data
#' @export
aggregate_samples <- function(TE,
                              metadata = c("Ribobase_QC_dedup_data"),
                              group_by = c("study", "cell_line"),
                              fun = c("mean")) {
  TE <- as.matrix(TE)
  if (is.null(colnames(TE))) {
    stop("TE must have column names matching Experiment identifiers")
  }

  group_choice <- match.arg(group_by)
  fun_choice <- match.arg(fun)

  utils::data(metadata, package = "ribobaser", envir = environment())
  metadata <- get(metadata, envir = environment())

  meta_subset <- metadata[!duplicated(metadata$Experiment),
                          c("Experiment", "Study", "Cell line"),
                          drop = FALSE]

  group_column <- if (group_choice == "study") "Study" else "Cell line"
  group_lookup <- stats::setNames(meta_subset[[group_column]], meta_subset$Experiment)
  sample_groups <- group_lookup[colnames(TE)]

  if (any(is.na(sample_groups))) {
    missing_samples <- colnames(TE)[is.na(sample_groups)]
    stop("Missing metadata for samples: ", paste(unique(missing_samples), collapse = ", "))
  }

  group_levels <- unique(sample_groups)

  if (identical(fun_choice, "mean")) {
    aggregated <- vapply(group_levels, function(g) {
      cols <- sample_groups == g
      rowMeans(TE[, cols, drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(TE)))
  } else {
    stop("Unsupported aggregation function: ", fun_choice)
  }

  dimnames(aggregated) <- list(rownames(TE), group_levels)
  aggregated
}
