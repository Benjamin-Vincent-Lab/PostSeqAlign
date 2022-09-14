human_ensembl_to_hgnc_entrez v0.1-01
https://sc.unc.edu/benjamin-vincent-lab/gene-sigs/human_ensembl_to_hgnc_entrez

This script was made to make a table that we can use to convert ENST transcripts 
to HGNC and entrez gene ids. cDNA transcript features have also been added.

##### Conversions
Using biomaRt and AnnotationDbi is a pretty awful way to go about converting 
genes. The resources are slow, lose connections, are down a lot, and have tons 
of gaps in the data. Contrarily, the GTF file has most of everything we need.  
To get this information we have to scrape it ourselves (GTF extraction tools 
won't get it :-/). If we combine that with the grch38 cDNA and we have 
everything to perfectly map the hgnc to the transcripts and pull out cDNA 
features. Entrez IDs are less complete. Entrez were better for matching to UCSC 
transcripts a while back, but aren't near as complate as HGNC anymore. At this
point it's highly recommended to map reads to ENST and convert to HGNC from 
there.

##### cDNA
In addition to doing the conversion of ensembl to different genes this script 
also gets the cDNA GC content and length.  The biomaRt `tx_length` and 
`percentage_gc_content` are for the full length transcript including introns.  
`cds_length` *should* be the length of the coding sequence (sum of exon 
lengths), but it's not even clsoe to the cDNA length.  The full length cDNA is 
not included in the final matrix to reduce the size of the table that will go 
into the R package but can be found in the cDNA intermdiate file.

##### Type and description
Information on the types (biomaRt `name_1006`) and description are provided for 
look-up separated by pipes so they can be easily searched for piped terms (eg 
'|AB|' will not show up in '|ABAB|')

##### UCSC caveats
USCS transcript ids, while on the table, were not the principal focus here so 
not a lot of effort was made to make those mappings 100% complete.  Plus adding 
UCSC would have duplicated a lot of ENST rows and made mapping to our principal 
goal less clear. 