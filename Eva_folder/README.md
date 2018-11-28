# Code files used in the manuscript  **xxxxx**      
___

These scripts form a pipeline for processing metabarcoding data from next-generation sequencing platforms. The pipeline consists of one bash script and three R scripts:  

 *  DADA2_demultiplex_tail_trim.sh  Demultiplexing and trimming of primers and tags.
 *  DADA2_start_181030.r   Basic quality filtering, chimera removal, and filtering of sequences likely to be caused by PCR and sequencing errors
 *  DADA2_end_181030.r   Combining data from different sequencing lanes and libraries, and outputting a list of OTUs and an OTU table of occurrence
 *  Tobias_taxonomy_v3_181030.r   Taxonomic classification of output from BLAST search of the OTU sequences. Uses NCBI GenBanks taxonomy and takes into account the top BLAST hits down to a user-specified similarity threshold."  
 
 Example files to run in with the scipts can be found in the directory above.  