
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' convert_ensembl
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ToDos
#' Make this take start from different types(eg: go from ensg to hgnc etc)
#' 
#' @title convert_ensembl 
#' 
#' @description 
#' Converts an Ensemble data.table into other gene types using a conversion table
#'
#' @param conversion_table_path Table to convert ensembl to other types
#' @param convert_to_column If using \code{get_human_ensembl_to_hgnc_entrez_path}, choose from: gene_id, hgnc_symbol, entrez_id or hgnc_entrez (eg, hgnc_symbol|entrez_id; not recommended, gene counts will be used several times when HGNC1|1 & HGNC2|1 map to the same ENST)
#' @param function_for_combining_counts Name of function to use to combine the genes. FWIW Joel Parker recommends that transcript colunt columns should be added to get gene level counts and microarrays should be averaged.
#' @param gene_biotypes The type of biomaRt gene_biotypes that should be output to the gene level output. Recommend you keep all of them. Many gene signature use lots of different biotypes.  Our bgvlab sigs use 18 types which is 228322/234393 of the ENST and 54367/59453 of the HGNC.
#' @param my_dt Input data.table with transcript counts. Should have samples by row, first columns should be the sample id with all subsequent columns headers in the Ensembl ENST format (no version). Assumes the first column is the sample_key.
#' @param readme_path Optional path to which the comments will be appended.
#' @param this_script_path Path to script that runs this function for documentation purposes
#' @return Returns a data.table with isoforms converted to genes.
#' 
#' @export
convert_ensembl = function(
  my_dt,
  conversion_table_path = get_human_ensembl_to_hgnc_entrez_path(),
  convert_to_column = "hgnc_symbol",
  gene_biotypes = NULL,
  function_for_combining_counts = sum,
  readme_path = NULL,
  this_script_path = NULL,
  thread_num = 1
){
  library(binfotron)
  library(magrittr)
  library(data.table)
  library(doMC)
  registerDoMC(thread_num)
  
  text_output = c(attributes(my_dt)$comments)
  
  if( convert_to_column == "hgnc_entrez" ){
    message("Highly recommend using individual entrez and hgnc matrices than to have both in one column.  We don't have 1:1 mappings for hgnc/entrez so this results in a lot of reads added to multiple times.")
  }
  
  if(!is.null(readme_path) && file.exists(readme_path)) file.remove(readme_path)
  
  a = function(...){
    my_output = paste0(...)
    assign('text_output', c(text_output, my_output), env = parent.frame())
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    cat(paste0(my_output,"\n"))
  }
  
  a("Running convert_ensembl")
  if(!is.null(this_script_path)) a("script location:", this_script_path)
  a(paste0("Converting transcript_id's (ie. ENST's) to ", convert_to_column, "'s."))
  
  a("")
  
  # convert to a data table according to the names of the list 
  #   (ie, use.names = T won't assume the genes are in the same order)
  sample_key = names(my_dt)[1]
  
  conversion_df = fread(conversion_table_path, data.table = F, stringsAsFactors = F)#, select = c("transcript_id","cdna_length"))
  conversion_df = conversion_df[!is.na(conversion_df[[convert_to_column]]),]
  conversion_df = conversion_df[conversion_df[[convert_to_column]] != "",]
  conversion_df = conversion_df[conversion_df$transcript_id %in% names(my_dt)[-1],]
  
  a(paste0("gene_biotypes: ", paste0(sort(unique(conversion_df$gene_biotype)), collapse = ", ")))
  a("")
  if(!is.null(gene_biotypes)){
    if(length(gene_biotypes) > 0){
      a(paste0("Restricting gene_biotypes to: ", paste0(gene_biotypes, collapse = ", ")))
      conversion_df = conversion_df[conversion_df$gene_biotype %in% gene_biotypes,]
    }
  }
  a("")
  start_trans_count = ncol(my_dt)-1
  my_dt = my_dt[, c(sample_key, names(my_dt)[names(my_dt) %in% conversion_df$transcript_id]), with =F] # went from 190K to 130K probably from loosing gene biotypes
  end_trans_count = ncol(my_dt)-1

  a(end_trans_count, "/", start_trans_count, " transcript id's were found in the conversion table.")

  my_genes = sort(unique(conversion_df[[convert_to_column]]))
  my_genes = my_genes[my_genes != ""]
  my_genes = my_genes[!is.na(my_genes)]
  # my_genes = my_genes[1:100000]
  # ptm <- proc.time()
  # my_dt = my_dt[,1, drop = FALSE]
  # for (my_gene in my_genes) {
  #   # get the transcripts it matches to
  #   my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
  #   my_dt[[my_gene]] = apply(my_dt[,my_transcripts, with = F],1,function_for_combining_counts)
  # }
  # message("for")
  # proc.time() - ptm
  # 
  # ptm <- proc.time()
  # list_counts = lapply(my_genes,function(my_gene){
  #   # get the transcripts it matches to
  #   my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
  #   return(apply(my_dt[,my_transcripts, with = F],1,function_for_combining_counts))
  # })
  # my_dt2 = cbind(my_dt[,1, drop = FALSE], data.frame(list_counts))
  # names(my_dt2)[2:ncol(my_dt2)] = my_genes
  # message("lapply")
  # proc.time() - ptm
  
  # ptm <- proc.time()
  # list_counts = parallel::mclapply(my_genes,function(my_gene){
  #   # get the transcripts it matches to
  #   my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
  #   return(apply(my_dt[,my_transcripts, with = F],1,function_for_combining_counts))
  # }, mc.cores = thread_num)
  # my_dt2 = cbind(my_dt[,1, drop = FALSE], data.frame(list_counts))
  # names(my_dt2)[2:ncol(my_dt2)] = my_genes
  # message("mclapply")
  # proc.time() - ptm
  
  #ptm <- proc.time()
  converted_dt = foreach(my_gene=my_genes, .combine=cbind) %dopar% { # 3x faster than for loop & lapply; 5% faster than mclapply
    my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
    apply(my_dt[,my_transcripts, with = F],1,function_for_combining_counts)
  }
  converted_dt = data.table(Run_ID = my_dt[[1]], converted_dt)
  
  names(converted_dt) = c(sample_key, my_genes)
  
  # message("foreach")
  #proc.time() - ptm
  
  attributes(converted_dt)$comments = c(text_output, "") 
  
  return(converted_dt)
}
