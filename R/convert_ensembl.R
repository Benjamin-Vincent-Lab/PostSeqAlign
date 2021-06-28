
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' convert_ensembl
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ToDos
#' Make this take start from different types(eg: go from ensg to hgnc etc)
#' 
#' @title convert_ensembl 
#' 
#' @description 
#' Converts an Ensemble matrix into other gene types using a conversion table

#' @param conversion_table_path Table to convert ensembl to other types
#' @param convert_to_column If using \code{get_human_ensembl_to_hgnc_entrez_path}, choose from: gene_id, hgnc_symbol, entrez_id or hgnc_entrez (eg, hgnc_symbol|entrez_id; not recommended, gene counts will be used several times when HGNC1|1 & HGNC2|1 map to the same ENST)
#' @param file_prefix Output file will have this prefix string appended to it.
#' @param function_for_combining_counts Name of function to use to combine the genes. FWIW Joel Parker recommends that transcript colunt columns should be added to get gene level counts and microarrays should be averaged.
#' @param gene_biotypes The type of biomaRt gene_biotypes that should be output to the gene level output. Recommend you keep all of them. Many gene signature use lots of different biotypes.  Our bgvlab sigs use 18 types which is 228322/234393 of the ENST and 54367/59453 of the HGNC.
#' @param input_file_path String path to input file. File should have samples by row, first columns should be the sample id with all subsequent columns headers in the Ensembl ENST format (no version). Assumes the first column is the sample_key.

#' @param output_dir Path to the output folder.
#' @param this_script_path Path to script that runs this function for documentation purposes
#' @return A vector of paths to the output file.
#' 
#' @export
convert_ensembl = function(
  conversion_table_path = get_human_ensembl_to_hgnc_entrez_path()
  convert_to_column = "hgnc_symbol",
  file_name = "converted_counts.tsv",
  gene_biotypes = NULL,
  input_file_path,
  method_for_combining_counts = sum,
  output_dir,
  this_script_path = '',
  thread_num = 1
){
  library(binfotron)
  library(magrittr)
  library(data.table)
  library(doMC)
  registerDoMC(thread_num)
  
  if(convert_to_column=="hgnc_entrez"){
    message("Highly recommend using individual entrez and hgnc matrices than to have both in one column.  We don't have 1:1 mappings for hgnc/entrez so this results in a lot of reads getting counted multiple times.")
  }
  
  dir.create(output_dir, showWarnings = F)
  
  readme_path = file.path(output_dir, paste0(file_prefix,"readme.txt"))
  if(file.exists(readme_path)){ file.remove(readme_path)}
  
  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    cat(paste0(my_output,"\n"))
  }
  
  a("")
  
  # convert to a data table according to the names of the list 
  #   (ie, use.names = T won't assume the genes are in the same order)
  dat = fread(input_file_path, stringsAsFactors = F)
  sample_key = names(dat)[1]
  
  conversion_df = fread(conversion_table_path, data.table = F, stringsAsFactors = F)#, select = c("transcript_id","cdna_length"))
  conversion_df = conversion_df[!is.na(conversion_df[[convert_to_column]]),]
  conversion_df = conversion_df[conversion_df[[convert_to_column]] != "",]
  
  a(paste0("* gene_biotypes: ", paste0(sort(unique(conversion_df$gene_biotype)), collapse = ", ")))
  if(!is.null(gene_biotypes)){
    if(length(gene_biotypes) > 0){
      a(paste0("* Restricting gene_biotypes to: ", paste0(gene_biotypes, collapse = ", ")))
      conversion_df = conversion_df[conversion_df$gene_biotype %in% gene_biotypes,]
    }
  }
  
  start_trans_count = ncol(dat)-1
  dat = dat[, c(sample_key, names(dat)[names(dat) %in% conversion_df$transcript_id]), with =F] # went from 190K to 130K probably from loosing gene biotypes
  end_trans_count = ncol(dat)-1

  a(end_trans_count, "/", start_trans_count, " transcript id's were found in the conversion table.")

  my_genes = sort(unique(conversion_df[[convert_to_column]]))
  my_genes = my_genes[my_genes != ""]
  my_genes = my_genes[!is.na(my_genes)]
  
  # my_genes = my_genes[1:100000]
  # ptm <- proc.time()
  # my_dt = dat[,1, drop = FALSE]
  # for (my_gene in my_genes) {
  #   # get the transcripts it matches to
  #   my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
  #   my_dt[[my_gene]] = apply(dat[,my_transcripts, with = F],1,function_for_combining_counts)
  # }
  # message("for")
  # proc.time() - ptm
  # 
  # ptm <- proc.time()
  # list_counts = lapply(my_genes,function(my_gene){
  #   # get the transcripts it matches to
  #   my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
  #   return(apply(dat[,my_transcripts, with = F],1,function_for_combining_counts))
  # })
  # my_dt2 = cbind(dat[,1, drop = FALSE], data.frame(list_counts))
  # names(my_dt2)[2:ncol(my_dt2)] = my_genes
  # message("lapply")
  # proc.time() - ptm
  
  # ptm <- proc.time()
  # list_counts = parallel::mclapply(my_genes,function(my_gene){
  #   # get the transcripts it matches to
  #   my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
  #   return(apply(dat[,my_transcripts, with = F],1,function_for_combining_counts))
  # }, mc.cores = thread_num)
  # my_dt2 = cbind(dat[,1, drop = FALSE], data.frame(list_counts))
  # names(my_dt2)[2:ncol(my_dt2)] = my_genes
  # message("mclapply")
  # proc.time() - ptm
  
  #ptm <- proc.time()
  my_dt = foreach(my_gene=my_genes, .combine=cbind) %dopar% { # 3x faster than for loop & lapply; 5% faster than mclapply
    my_transcripts = unique(conversion_df$transcript_id[conversion_df[[convert_to_column]] == my_gene])
    apply(dat[,my_transcripts, with = F],1,function_for_combining_counts)
  }
  my_dt = cbind(dat[,1, drop = FALSE], my_dt)
  names(my_dt)[2:ncol(my_dt)] = my_genes
  # message("foreach")
  #proc.time() - ptm
  
  
  output_path = file.path(output_dir, paste0(file_name))
  
  fwrite(dat, output_path, sep = "\t")

}
