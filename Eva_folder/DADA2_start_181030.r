library(dada2)
setwd("your_wd")
pwd <- getwd()
setwd(pwd)
main_path <- pwd

#This script is run for each library in each sequencing lane. Put the resulting seqtab files in 
#separate folders for each library, naming the files to reflect the sequencing lane, e.g. 
#"seqtab2_AS" for the file from lane 2. Then continue with the script DADA2_end.r

#Matching the "Sense"-reads.
path <- file.path(main_path, "DADA2_SS")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("-R1.", fastqs)]
fnRs <- fastqs[grepl("-R2.", fastqs)]
sample.names <- sapply(strsplit(fnFs, "-"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "matched")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
for(i in seq_along(fnFs)) {
      if (file.info(fnFs[i])$size == 0) {
      print(paste(fnFs[i], "empty.")) } else {
        print(paste("processing", fnFs[i]))
        fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                          matchIDs=TRUE)
      }
  }

#Matching the "Anti-Sense"-reads.
path <- file.path(main_path, "DADA2_AS")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("-R2.", fastqs)] # Reverse direction compared to above for the "sense reads". In practice the reads are here complement reversed to be in the same orientation as the "sense" reads.
fnRs <- fastqs[grepl("-R1.", fastqs)] # See above
sample.names <- sapply(strsplit(fnFs, "-"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "matched")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
for(i in seq_along(fnFs)) {
      if (file.info(fnFs[i])$size == 0) {
      print(paste(fnFs[i], "empty.")) } else {
        print(paste("processing", fnFs[i]))
        fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                          matchIDs=TRUE)
      }
    }

#filtering of the Sense-reads:
path <- file.path(main_path, "DADA2_SS/matched") 
fns <- list.files(path)
fastqs <- fns[grepl("matched.fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
for(i in seq_along(fnFs)) {
  if (file.info(fnFs[i])$size == 0) {
    print(paste(fnFs[i], "empty.")) } else {
      print(paste("processing", fnFs[i]))
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    minLen=10, maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
  }
}
#filtering of the antiSense-reads:
path <- file.path(main_path, "DADA2_AS/matched") 
fns <- list.files(path)
fastqs <- fns[grepl("matched.fastq", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]), 
                    minLen=10, maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE)
}

#Processing the set of files containing the forward primer in the R1 reads (the sense reads):
filt_path <- file.path(main_path, "DADA2_SS/matched/filtered") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE)
dadaRs <- dada(derepRs, err=NULL, selfConsist = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_SS <- makeSequenceTable(mergers[names(mergers)])
seqtab.nochim_SS <- removeBimeraDenovo(seqtab_SS, verbose=TRUE)
stSS <- file.path(main_path,"seqtab_SS")
stnsSS <- file.path(main_path,"seqtab.nochim_SS")
saveRDS(seqtab_SS,stSS)
saveRDS(seqtab.nochim_SS,stnsSS)

#Then DADA2 processing of "the antisense" reads:
filt_path <- file.path(main_path, "DADA2_AS/matched/filtered") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE)
dadaRs <- dada(derepRs, err=NULL, selfConsist = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_AS <- makeSequenceTable(mergers[names(mergers)])
seqtab.nochim_AS <- removeBimeraDenovo(seqtab_AS, verbose=TRUE)
stAS <- file.path(main_path,"seqtab_AS")
stnsAS <- file.path(main_path,"seqtab.nochim_AS")
saveRDS(seqtab_AS,stAS)
saveRDS(seqtab.nochim_AS,stnsAS)

