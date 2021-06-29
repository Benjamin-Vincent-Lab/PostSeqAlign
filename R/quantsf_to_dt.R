
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' quantsf_to_dt
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title quantsf_to_dt 
#' 
#' @description 
#' Joins individual sample files into one tsv file.
#' 
#' @param counts_or_tpm 'counts' or 'tpm'
#' @param input_file_paths Named character vector of paths to the pipeline output data. Samples will be named by the names of the vector. If not named the will be named after the quant file.
#' @param readme_path Optional path to which the comments will be appended.
#' @param sample_key String to name the sample identifier column.
#' @param this_script_path Path to script that runs this function for documentation purposes
#' @return A vector of paths to the output file.
#' 
#' @export
quantsf_to_dt = function(
  input_file_paths,# = system(paste0("ls ", RAW_DATA_DIR, "/pipeline_output/star_salmon/*/*_quant.sf"), intern = TRUE)
  counts_or_tpm = "counts",
  readme_path = NULL,
  sample_key = "Run_ID",
  this_script_path = ''
){
  library(magrittr)
  library(data.table)
  
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
  
  a("Running quantsf_to_dt")
  if(!is.null(this_script_path)) a("script: ", this_script_path)
  a("")
  a(paste0("Assembling quant.sf ", counts_or_tpm, " into data.table."))
  a("")
  
  if(length(names(input_file_paths)) == 0){
    names(input_file_paths) = gsub(".quant.sf$", "", basename(input_file_paths))
  }
  
  
  counts_or_tpm = tolower(counts_or_tpm)
  if(counts_or_tpm == "counts"){
    salmon_column = "NumReads"
  } else if (counts_or_tpm == "tpm"){
    salmon_column = "TPM"
  }
  
  a("Reading in files from input_file_paths:")
  for(input_file_path in input_file_paths){a("  ", input_file_path)}
  read_data = lapply(input_file_paths, function(input_file_path){
    counts_df = fread(input_file_path, select = c("Name", salmon_column), data.table = F)
    return_list = lapply(counts_df[[salmon_column]], function(x)x) # rbindlist needs a list so we turn this into a list
    names(return_list) = counts_df[["Name"]] # Name the list items so they get assigned to the right column
    return(return_list)
  })
  a("")
  
  # convert to a data table according to the names of the list 
  #   (ie, use.names = T won't assume the genes are in the same order)
  dat = rbindlist(read_data, use.names = TRUE, fill = TRUE)
  
  # add the sample names
  dat = data.table(Run_ID = names(read_data), dat)
  if (sample_key != "Run_ID") names(dat)[names(dat) == "Run_ID"] = sample_key
  
  row.names(dat) = NULL
  
  # some grch38 references (eg GATK) put decimals at the end of the names. 
  names(dat)[2:ncol(dat)] %<>% substring(., 1, 15) 
  attributes(dat)$comments = c(text_output, "") 
  
  return(dat)
}
