library(dada2)

#Define a function for combining two or more tables, collapsing samples with similar names:  
sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
  # Combine passed tables into a list
  tables <- list(table1, table2)
  tables <- c(tables, list(...))
  # Validate tables
  if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
    stop("At least two valid sequence tables, and no invalid objects, are expected.")
  }
  sample.names <- rownames(tables[[1]])
  for(i in seq(2, length(tables))) {
    sample.names <- c(sample.names, rownames(tables[[i]]))
  }
  seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
  sams <- unique(sample.names)
  # Make merged table
  rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
  rownames(rval) <- sams
  colnames(rval) <- seqs
  for(tab in tables) {
    rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
  }
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
  }
  rval
}

setwd("your_wd_library_1")
#setwd("your_wd_library_2")
#setwd("your_wd_library_3")
#setwd("your_wd_library_4")
pwd <- getwd()
setwd(pwd)
main_path <- pwd

###The below (until line 106) is done separately for each library, changing the working directory as indicated above
###Lane 1
stAS <- file.path(main_path,"seqtab_AS")
stnsAS <- file.path(main_path,"seqtab.nochim_AS")
stSS <- file.path(main_path,"seqtab_SS")
stnsSS <- file.path(main_path,"seqtab.nochim_SS")
seqtab.nochim_AS <- readRDS(stnsAS)
seqtab.nochim_SS <- readRDS(stnsSS)
seqtab_AS <- readRDS(stAS)
seqtab_SS <- readRDS(stSS)
sumtable <- sumSequenceTables(seqtab_SS,seqtab_AS) # den her kommando kan merge to tabeller. sammen navn og/eller samme sekvens merges.
nochim_sumtable <- sumSequenceTables(seqtab.nochim_SS,seqtab.nochim_AS)

###Lane 2
stAS2 <- file.path(main_path,"seqtab2_AS")
stnsAS2 <- file.path(main_path,"seqtab2.nochim_AS")
stSS2 <- file.path(main_path,"seqtab2_SS")
stnsSS2 <- file.path(main_path,"seqtab2.nochim_SS")
seqtab.nochim_AS2 <- readRDS(stnsAS2)
seqtab.nochim_SS2 <- readRDS(stnsSS2)
seqtab_AS2 <- readRDS(stAS2)
seqtab_SS2 <- readRDS(stSS2)
sumtable2 <- sumSequenceTables(seqtab_SS2,seqtab_AS2) # den her kommando kan merge to tabeller. sammen navn og/eller samme sekvens merges.
nochim_sumtable2 <- sumSequenceTables(seqtab.nochim_SS2,seqtab.nochim_AS2)

###Lane 3
stAS3 <- file.path(main_path,"seqtab3_AS")
stnsAS3 <- file.path(main_path,"seqtab3.nochim_AS")
stSS3 <- file.path(main_path,"seqtab3_SS")
stnsSS3 <- file.path(main_path,"seqtab3.nochim_SS")
seqtab.nochim_AS3 <- readRDS(stnsAS3)
seqtab.nochim_SS3 <- readRDS(stnsSS3)
seqtab_AS3 <- readRDS(stAS3)
seqtab_SS3 <- readRDS(stSS3)
sumtable3 <- sumSequenceTables(seqtab_SS3,seqtab_AS3) # den her kommando kan merge to tabeller. sammen navn og/eller samme sekvens merges.
nochim_sumtable3 <- sumSequenceTables(seqtab.nochim_SS3,seqtab.nochim_AS3)

###Lane 4
stAS4 <- file.path(main_path,"seqtab4_AS")
stnsAS4 <- file.path(main_path,"seqtab4.nochim_AS")
stSS4 <- file.path(main_path,"seqtab4_SS")
stnsSS4 <- file.path(main_path,"seqtab4.nochim_SS")
seqtab.nochim_AS4 <- readRDS(stnsAS4)
seqtab.nochim_SS4 <- readRDS(stnsSS4)
seqtab_AS4 <- readRDS(stAS4)
seqtab_SS4 <- readRDS(stSS4)
sumtable4 <- sumSequenceTables(seqtab_SS4,seqtab_AS4) # den her kommando kan merge to tabeller. sammen navn og/eller samme sekvens merges.
nochim_sumtable4 <- sumSequenceTables(seqtab.nochim_SS4,seqtab.nochim_AS4)

###Merge Lanes
sumtable_1_2 <- sumSequenceTables(sumtable,sumtable2) # den her kommando kan merge to tabeller. sammen navn og/eller samme sekvens merges.
nochim_sumtable_1_2 <- sumSequenceTables(nochim_sumtable,nochim_sumtable2)
sumtable_3_4 <- sumSequenceTables(sumtable3,sumtable4) # den her kommando kan merge to tabeller. sammen navn og/eller samme sekvens merges.
nochim_sumtable_3_4 <- sumSequenceTables(nochim_sumtable3,nochim_sumtable4)
sumtable_1_4 <- sumSequenceTables(sumtable_1_2,sumtable_3_4) # den her kommando kan merge to tabeller. sammen navn og/eller samme sekvens merges.
nochim_sumtable_1_4 <- sumSequenceTables(nochim_sumtable_1_2,nochim_sumtable_3_4)

stBoth <- file.path(main_path,"seqtab_Both")
stnsBoth <- file.path(main_path,"seqtab.nochim_Both")
saveRDS(sumtable_1_4,stBoth)
saveRDS(nochim_sumtable_1_4,stnsBoth)

###In the below, the libraries are merged
setwd("your_wd_library_1")
table_library_1 <- readRDS("seqtab.nochim_Both")
setwd("your_wd_library_2")
table_library_2 <- readRDS("seqtab.nochim_Both")
setwd("your_wd_library_3")
table_library_3 <- readRDS("seqtab.nochim_Both")
setwd("your_wd_library_4")
table_library_4 <- readRDS("seqtab.nochim_Both")

sum_1_2 <- sumSequenceTables(table_library_1,table_library_2)
sum_1_2_3 <- sumSequenceTables(sum_1_2,table_library_3)
sum_1_2_3_4 <- sumSequenceTables(sum_1_2_3,table_library_4)

#Transpose table, assign names, extract sequences and saving table, for further processing:
trasnochim_sumtable <- as.data.frame(t(sum_1_2_3_4))
#Get DNA sequences
sequences <- row.names(trasnochim_sumtable)
#Assign new rownames
row.names(trasnochim_sumtable) <- paste0("seq",seq.int(nrow(trasnochim_sumtable)))
tbname <- file.path(main_path,"DADA2_nochim.table")
{write.table(trasnochim_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}
#Extract OTUs (sequences)
sinkname <- file.path(main_path,"DADA2_nochim.otus")
{
  sink(sinkname)
  for (seqX in seq.int(nrow(trasnochim_sumtable))) {
    header <- paste0(">","seq",seqX,"\n")
    cat(header)
    seqq <- paste0(sequences[seqX],"\n")
    cat(seqq)
  }
  sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table) {
  out_path <- file.path(main_path, "DADA2_extracted_samples_nochim")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",rownames(my_table)[sampleX],";size=",my_table[sampleX,seqX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(sum_1_2_3_4)

#The process above is repeated for the raw sequences. 
setwd("your_wd_library_1")
table_library_1 <- readRDS("seqtab_Both")
setwd("your_wd_library_2")
table_library_2 <- readRDS("seqtab_Both")
setwd("your_wd_library_3")
table_library_3 <- readRDS("seqtab_Both")
setwd("your_wd_library_4")
table_library_4 <- readRDS("seqtab_Both")

sum_1_2 <- sumSequenceTables(table_library_1,table_library_2)
sum_1_2_3 <- sumSequenceTables(sum_1_2,table_library_3)
sum_1_2_3_4 <- sumSequenceTables(sum_1_2_3,table_library_4)

#Transpose table, assign names, extract sequences and saving table, for further processing:
trasraw_sumtable <- as.data.frame(t(sum_1_2_3_4))
#Get DNA sequences
sequences <- row.names(trasraw_sumtable)
#Assign new rownames
row.names(trasraw_sumtable) <- paste0("seq",seq.int(nrow(trasraw_sumtable)))
tbname <- file.path(main_path,"DADA2_raw.table")
{write.table(trasraw_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}
#Extract OTUs (sequences)
sinkname <- file.path(main_path,"DADA2_raw.otus")
{
  sink(sinkname)
  for (seqX in seq.int(nrow(trasraw_sumtable))) {
    header <- paste0(">","seq",seqX,"\n")
    cat(header)
    seqq <- paste0(sequences[seqX],"\n")
    cat(seqq)
  }
  sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table) {
  out_path <- file.path(main_path, "DADA2_extracted_samples_raw")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",rownames(my_table)[sampleX],";size=",my_table[sampleX,seqX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(sum_1_2_3_4)
