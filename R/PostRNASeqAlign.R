#' PostRNASeqAlign: Functions for writing commands and post processing the Star Salmon pipeline.
#'
#' Given a sample data matrix that indicate the sample name, sample folder name and input paths
#' this script can write out pipeline commands to run star salmon on the cluster and post process the
#' data into ucsc_isoform_counts.tsv and hgnc_entrezid_gene_counts.tsv. 
#' 
#' @import data.table parallel annotate binfotron
#' 
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"