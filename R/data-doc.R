#' Ribo count matrix
#'
#' Genes by samples count matrix.
#' @format A numeric matrix with G genes and S samples
#' @source Yue's Mozart Folder. Assuming that this follows Jonathan's counting
"ribo_raw_human_cap_995"

#' RNA count matrix
#' @format A numeric matrix
#' @source Yue's Mozart Folder. Assuming that this follows Jonathan's counting
"rnaseq_raw_human_cap_995"

#' RiboBase dedup counts QC metrics
#' @format A data.frame with 14 variables
#' \describe{
#'   \item{Experiment}{Unique identifier for the experiment GSM (character).}
#'   \item{Study}{Indicates Study ID GSE (character).}
#'   \item{Cell line}{Cell type (character).}
#'   \item{Species}{human/mouse (character).}
#'   \item{Start Length}{Dynamic Cutoff RPF length start (num).}
#'   \item{End Length}{Dynamic Cutoff RPF length end (num).}
#'   \item{Periodicitiy Score}{The weighted sum of dominant frame percentages, with weights based on the total read counts at each length(num).}
#'   \item{Periodicity distr}{The percentage of reads mapped to the dominant frame across the dynamic RPF range(chr).}
#'   \item{CDS coverage 15_40}{Percent Coverage (num).}
#'   \item{CDS coverage 27_30}{Percent Coverage (num).}
#'   \item{CDS coverage dynamic}{Percent Coveraget (num).}
#'   \item{Read counts 15_40}{CDS mapping reads in range (num).}
#'   \item{Read counts 27_30}{CDS mapping reads in range (num).}
#'   \item{Read counts dynamic}{CDS mapping reads in range (num).}
#' }
#' @source Hurley.
"Ribobase_QC_dedup_data"

#' RiboBase dedup counts QC metrics
#' @format A data.frame with 14 variables
#' \describe{
#'   \item{Experiment}{Unique identifier for the experiment GSM (character).}
#'   \item{Study}{Indicates Study ID GSE (character).}
#'   \item{Cell line}{Cell type (character).}
#'   \item{Species}{human/mouse (character).}
#'   \item{Start Length}{Dynamic Cutoff RPF length start (num).}
#'   \item{End Length}{Dynamic Cutoff RPF length end (num).}
#'   \item{Periodicitiy Score}{The weighted sum of dominant frame percentages, with weights based on the total read counts at each length(num).}
#'   \item{Periodicity distr}{The percentage of reads mapped to the dominant frame across the dynamic RPF range(chr).}
#'   \item{CDS coverage 15_40}{Percent Coverage (num).}
#'   \item{CDS coverage 27_30}{Percent Coverage (num).}
#'   \item{CDS coverage dynamic}{Percent Coveraget (num).}
#'   \item{Read counts 15_40}{CDS mapping reads in range (num).}
#'   \item{Read counts 27_30}{CDS mapping reads in range (num).}
#'   \item{Read counts dynamic}{CDS mapping reads in range (num).}
#' }
#' @source Hurley.
"Ribobase_QC_non_dedup_data"
