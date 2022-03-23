
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ensembl_counts_to_rkpm
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title ensembl_counts_to_rkpm 
#' 
#' @description 
#' Converts an Ensemble counts matrix into RKPM. <rant> RKPM is a terrible way to 
#' look at your data.  Only do it if you have to. If you choose cDNA length in 
#' your calculation then you loose all the non-coding transcripts.  If you use 
#' the full transcript length then you should account for the number of full and 
#' cds sequences which there's no way to know...</rant>
#' @param my_dt data.table with samples by row, first column should be the sample id with all subsequent columns headers in the Ensembl ENST format (no version).
#' @param lengths_table_path Path to table for with cdna_length and transcript_id columns.
#' @param readme_path Optional path to which the comments will be appended.
#' @param this_script_path Path to script that runs this function for documentation purposes
#' @return A vector of paths to the output file.
#' 
#' @export
ensembl_counts_to_rkpm = function(
  my_dt,
  lengths_table_path = get_human_ensembl_to_hgnc_entrez_path(),
  readme_path = NULL,
  this_script_path = ''
){
  library(magrittr)
  library(data.table)
  library(matrixStats)
  
  text_output = c()
  
  if(!is.null(readme_path) && file.exists(readme_path)) file.remove(readme_path)
  
  a = function(...){
    my_output = paste0(...)
    assign('text_output', c(text_output, my_output), env = parent.frame())
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    cat(paste0(my_output,"\n"))
  }
  
  a("Running ensembl_counts_to_rkpm")
  if(!is.null(this_script_path)) a("script: ", this_script_path)
  a("")
  
  # convert to a data table according to the names of the list 
  #   (ie, use.names = T won't assume the genes are in the same order)
  sample_key = names(my_dt)[1]
  
  lengths_df = fread(lengths_table_path, data.table = F, select = c("transcript_id","cdna_length"))
  
  a("Droping transcripts where cDNA length is NA")
  lengths_df = lengths_df[!is.na(lengths_df$cdna_length), ]
  lengths_df = lengths_df[lengths_df$cdna_length > 0, ]
  
  a("Calculating cdna_length/1000.")
  transcript_length_per_k_lut = lengths_df$cdna_length/1000
  names(transcript_length_per_k_lut) = lengths_df$transcript_id
  rm(lengths_df)
  transcript_length_per_k_lut %<>% .[names(.) %in% names(my_dt)]
  my_dt %<>% .[,c(sample_key, names(transcript_length_per_k_lut)), with = F]
  
  # tried geting sample totals with data.frame operation but it was slower
  # sub_dat = dat[, 1:40000]
  # with 4 cores 40K columns ran in 19 sec instead of 26
  # with 24 cores all isoforms took 1130 sec
  # ptm <- proc.time()
  # sample_totals = dat[, parallel::mclapply(.SD, sum, na.rm=TRUE, mc.cores = thread_num), .SDcols=2:ncol(dat) ] 
  # proc.time() - ptm
  # use () around it to have data.table check the var name
  
  a("Calculating sample totals.")
  my_dt[, `:=`(total = rowSums(.SD, na.rm=T)), .SDcols=2:ncol(my_dt), by=1:nrow(my_dt)] # 56-60 sec

  total_lut = my_dt$total/1000000
  names(total_lut) = my_dt[[1]]
  my_dt[ , total := NULL]

  # https://stackoverflow.com/questions/20596433/how-to-divide-each-row-of-a-matrix-by-elements-of-a-vector-in-r
  mat = my_dt[,2:ncol(my_dt)] %>% as.matrix
  #ptm <- proc.time()
  a("Divide each gene column by its cdna_length/1000")
  mat = t(t(mat) / transcript_length_per_k_lut)
  
  a("Divide each sample row by its sample_total/1,000,000")
  mat = mat/total_lut
  
  #proc.time() - ptm
  my_dt[,2:ncol(my_dt)] = data.frame(mat)
  
  
  attributes(my_dt)$comments = c(text_output, "") 
  
  return(my_dt)

}
