
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' quantsf_to_matrix
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title quantsf_to_matrix 
#' 
#' @description 
#' Joins individual sample files into one tsv file.
#' 
#' @param counts_or_tpm 'counts' or 'tpm'
#' @param file_prefix Output files will have this prefix string appended to them
#' @param input_file_paths Named character vector of paths to the pipeline output data. Samples will be named by the names of the vector. If not named the will be named after the quant file.
#' @param output_dir Path to the output folder.
#' @param sample_key String to name the sample identifier column.
#' @param this_script_path Path to script that runs this function for documentation purposes
#' @return A vector of paths to the output file.
#' 
#' @export
quantsf_to_matrix = function(
  counts_or_tpm = "counts",
  file_prefix = "",
  input_file_paths,# = system(paste0("ls ", RAW_DATA_DIR, "/pipeline_output/star_salmon/*/*_quant.sf"), intern = TRUE)
  output_dir,# = file.path(base_dir, "post_processing", "star_salmon")
  sample_key = "Run_ID",
  this_script_path = ''
){
  library(binfotron)
  library(magrittr)
  library(data.table)
  
  dir.create(output_dir, showWarnings = F)
  
  if(file_prefix != ""){file_prefix %<>% paste0("__")}
  
  readme_path = file.path(output_dir, paste0(file_prefix,"readme.txt"))
  if(file.exists(readme_path)){ file.remove(readme_path)}
  
  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    cat(paste0(my_output,"\n"))
  }
  
  a(paste0("## Making isoform counts matrix: ", this_script_path))
  a("")
  
  if(length(names(input_file_paths)) == 0){
    names(input_file_paths) = gsub(".quant.sf$", "", basename(input_file_paths))
  }
  
  
  a("Reading in files from input_file_paths:")
  counts_or_tpm = tolower(counts_or_tpm)
  if(counts_or_tpm == "counts"){
    salmon_column = "NumReads"
  } else if (counts_or_tpm == "tpm"){
    salmon_column = "TPM"
  }
  
  read_data = lapply(input_file_paths, function(input_file_path){
    a("  ", input_file_path)
    counts_df = fread(input_file_path, select = c("Name", salmon_column), data.table = F)
    return_list = lapply(counts_df[[salmon_column]], function(x)x) # rbindlist needs a list so we turn this into a list
    names(return_list) = counts_df[["Name"]] # Name the list items so they get assigned to the right column
    return(return_list)
  })
  a("")
  
  # convert to a data table according to the names of the list 
  #   (ie, use.names = T won't assume the genes are in the same order)
  dat = rbindlist(read_data, use.names = TRUE, fill = TRUE)
  
  
  output_paths = c()
  
  # add the sample names
  dat = data.frame(Run_ID = names(read_data), dat)
  if (sample_key != "Run_ID") names(dat)[names(dat) == "Run_ID"] = sample_key
  
  row.names(dat) = NULL
  
  # some grch38 references (eg GATK) put decimals at the end of the names. 
  names(dat)[2:ncol(dat)] %<>% substring(., 1, 15) 
  
  isoform_output_path = file.path(output_dir, paste0(file_prefix, "trans_", counts_or_tpm,".tsv"))
  
  fwrite(dat, isoform_output_path, sep = "\t")
  output_paths = c(output_paths, isoform_output_path)
  
  return(output_paths)
}
