
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ensembl_counts_to_tpm
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title ensembl_counts_to_tpm 
#' 
#' @description 
#' Converts an Ensemble ENST counts matrix into TPM.
#' @param my_dt data.table with samples by row, first column should be the sample id with all subsequent columns headers in the Ensembl ENST format (no version).
#' @param lengths_table_path Path to table for with cdna_length and transcript_id columns.
#' @param readme_path Optional path to which the comments will be appended.
#' @param this_script_path Path to script that runs this function for documentation purposes
#' @return A vector of paths to the output file.
#' 
#' @export
ensembl_counts_to_tpm = function(
  my_dt,
  lengths_table_path = get_human_ensembl_to_hgnc_entrez_path(),
  readme_path = NULL,
  this_script_path = ''
){
  library(magrittr)
  library(data.table)
  library(matrixStats)
  
	# to manually inspect calculations, although this should have a unit test made for it
	# setwd("/datastore/nextgenout5/share/labs/Vincent_Lab/tools/raft/projects/VL340_SITC_2020/20220701_tissue_specific_disc_val/batch/transcript_counts/merged_enst_counts")
	# my_dt = data.table::fread("merged_enst_counts.tsv", select = 1:50)
	# my_dt = my_dt[grepl('^Prin',my_dt$Run_ID), ]
	# lengths_table_path = PostRNASeqAlign::get_human_ensembl_to_hgnc_entrez_path()
	
	# source https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
	# TPM is very similar to RPKM and FPKM. The only difference is the order of operations. Here’s how you calculate TPM:
	# Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
	# Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
	# Divide the RPK values by the “per million” scaling factor. This gives you TPM.
	
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
  transcript_length_in_kb_lut = lengths_df$cdna_length/1000
  names(transcript_length_in_kb_lut) = lengths_df$transcript_id
  rm(lengths_df)
  transcript_length_in_kb_lut %<>% .[names(.) %in% names(my_dt)]
  my_dt %<>% .[,c(sample_key, names(transcript_length_in_kb_lut)), with = F]
  
  # tried geting sample totals with data.frame operation but it was slower
  # sub_dat = dat[, 1:40000]
  # with 4 cores 40K columns ran in 19 sec instead of 26
  # with 24 cores all isoforms took 1130 sec
  # ptm <- proc.time()
  # sample_totals = dat[, parallel::mclapply(.SD, sum, na.rm=TRUE, mc.cores = thread_num), .SDcols=2:ncol(dat) ] 
  # proc.time() - ptm
  # use () around it to have data.table check the var name
 
  # https://stackoverflow.com/questions/20596433/how-to-divide-each-row-of-a-matrix-by-elements-of-a-vector-in-r
  #ptm <- proc.time()
  a("Divide each gene column by its cdna_length/1000")
  mat = my_dt[,2:ncol(my_dt)] %>% as.matrix
  mat = t(t(mat) / transcript_length_in_kb_lut)
  # manually confirmed by looking up transcript_length, dividing by 1000
  # dividing my_dt column, and comparing to mat
  my_dt[,2:ncol(my_dt)] = data.frame(mat)
  
  
  a("Calculating sample RPK totals.")
  my_dt[, `:=`(total = rowSums(.SD, na.rm=T)), .SDcols=2:ncol(my_dt), by=1:nrow(my_dt)] # 56-60 sec

  per_million_scaling_factor_lut = my_dt$total/1000000
  names(per_million_scaling_factor_lut) = my_dt[[1]]
  # > apply(my_dt[1:10,2:26], 1, sum)
  # [1]  80.68529  42.08124 188.85702  59.57138  66.66333 124.75849 200.12822  56.76575 141.43841 379.95485
  # > my_dt$total[1:10]
  # [1]  80.68529  42.08124 188.85702  59.57138  66.66333 124.75849 200.12822  56.76575 141.43841 379.95485
  my_dt[ , total := NULL]
  
  mat = my_dt[,2:ncol(my_dt)] %>% as.matrix
  
  a("Divide each sample row by its sample_total/1,000,000")
  mat = mat/per_million_scaling_factor_lut
  
  #proc.time() - ptm
  
  my_dt[,2:ncol(my_dt)] = data.frame(mat)
  
  attributes(my_dt)$comments = c(text_output, "") 
  
  return(my_dt)

}
