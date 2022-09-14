Post processes output from RNA-Seq alignments.  Has been used on salmon 'quant.sf' outputs, 
but may generalize to other outputs soon.



## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_github("Benjamin-Vincent-Lab/PostRNASeqAlign")
```

Or for a specific version:
``` r
devtools::install_github("Benjamin-Vincent-Lab/PostRNASeqAlign", ref = "0.5.0")
```

## Previous locations
https://github.com/Benjamin-Vincent-Lab/StarSalmon
https://sc.unc.edu/dbortone/starsalmon
https://sc.unc.edu/benjamin-vincent-lab/starsalmon
Moved to github so that the package could be accessed without a token.

## Replacing `post_process_rnaseq_align`
`post_process_rnaseq_align` was a patchwork of mapping isoforms to genes using 
biomaRt/AnnotationDbi and it still wasn't able to map all of the HGNC symbols 
to ENSTs one-to-one. A new conversion table has been added (see repo 
human_ensembl_to_hgnc_entrez v0.1-01). These new functions have been created to 
take advantage of this table:
* `get_human_ensembl_to_hgnc_entrez_path` - for accessing the conversion table
* `quantsf_to_dt` - just takes counts or tpm from the quant.sf and puts them in a matrix
* `convert_ensembl` - converts ENST's to ENSG, HGNC symobls or Entrez IDs
* `ensembl_counts_to_rkpm` - converts Ensemble counts to RKPM

These smaller functions no-longer upper-quartile normalization or log2 transform 
the data, but these transformations can be done using binfotron 
(`log_transform_plus` & `normalize_rows_by_quartile`).

Finally, using 8 threads, `convert_ensembl` is about 3x faster than 
`post_process_rnaseq_align`.

`post_process_rnaseq_align` will be kept around until it's phased out of our 
workflows.

## Recommendations
* Avoid using RKPM unless you are doing it preparation for other module (eg TIDE).
* It's better to stick with HGNC instead of Entrez IDs as HGNCs map one-to-one 
and the Entrez data have a lot of gaps.

##ToDos
* Right now the `human_ensembl_to_hgnc_entrez` table was built using a GTF for 
GRCh38/v103.  If a different GTF is used to map the reads then ideally a new 
conversion table would be provided for those conversions.
* Make unit tests
* Conversion tables should only convert 1 to 1 to avoid problems where one 
column's content could be duplicated if a one-to-many relationship exists in 
another column

## Example code
``` R
transcript_counts_dt = quantsf_to_dt(
  input_file_paths = list.files("~/_tagged_batches/Gide_Cell_2019/rna_quant/all_v1", pattern = "*quant.sf", full.names = T)[1:5],
  counts_or_tpm = "counts", # or tpm
  readme_path = "~/scratch/test_readme.txt",
  sample_key = "Run_ID",
  this_script_path = housekeeping::get_script_dir_path(include_file_name = T)
)

hgnc_counts_dt = convert_ensembl(
  transcript_counts_dt,
  conversion_table_path = PostRNASeqAlign::get_human_ensembl_to_hgnc_entrez_path(),
  convert_to_column = "hgnc_symbol", # gene_id, hgnc_symbol, entrez_id or hgnc_entrez
  gene_biotypes = NULL, # eg: protein_coding
  function_for_combining_counts = sum,
  readme_path = "~/scratch/test_readme2.txt",
  this_script_path = housekeeping::get_script_dir_path(include_file_name = T),
  thread_num = 8
)

rkpm_dt = ensembl_counts_to_rkpm(
  transcript_counts_dt,
  lengths_table_path = PostRNASeqAlign::get_human_ensembl_to_hgnc_entrez_path(),
  readme_path = "~/scratch/test_readme3.txt",
  this_script_path = housekeeping::get_script_dir_path(include_file_name = T)
)
```