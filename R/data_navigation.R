#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_biomart_hg38_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_biomart_hg38_path 
#' 
#' @description 
#' Returns the path a saved mart
#' 
#' @param none
#' 
#' @return A path to the rds file.
#' 
#' @family mart
#' 
#' @export
get_biomart_hg38_path = function(){
  return(system.file("biomart", "hg38", "hsa_ensembl_ucsc.tsv", package = "PostRNASeqAlign"))
}


#' get_human_ensembl_to_hgnc_entrez_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_human_ensembl_to_hgnc_entrez_path 
#' 
#' @param none
#' 
#' @return A path to the table for converting Ensembl ENSTs to HGNC and Entrez IDs
#' 
#' @export
get_human_ensembl_to_hgnc_entrez_path = function(){
  return(system.file("human_ensembl_to_hgnc_entrez", "transcript_conversion_table.tsv", package = "PostRNASeqAlign"))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_biomart_grch38_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_biomart_grch38_path 
#' 
#' @description 
#' Returns the path to a saved mart
#' 
#' @param none
#' 
#' @return A path to the rds file.
#' 
#' @family mart
#' 
#' @export
get_biomart_grch38_path = function(){
  return(system.file("biomart", "grch38", "bm_result.tsv", package = "PostRNASeqAlign"))
}

