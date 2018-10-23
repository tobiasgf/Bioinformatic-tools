#Functions for classifying OTUs
#
# Author: Tobias G. Frøslev 27.1.2018
# Modified substantioally 22.10.2018 - info/instructions not updated!

# requires a blast output from a set of OTUs
# minimum number of fields requires are qseqid, staxid, pident, ssciname and evalue
# optimal blastn command to use:
#
# blastn -remote -db nt -max_target_seqs 40 -outfmt "6 std qlen qcovs sgi sseq ssciname staxid" -out BLAST_HIT_OUTPUT -qcov_hsp_perc 90 -perc_identity 80 -query INPUT

# Instructions
# 1) Load the three functions below: assign_taxonomy, prefilter, get_classification, evaluate_classification
# 2) Classify you OTUs by running the wrapper function assign_taxonomy like this:
#   classified_table <- assign_taxonomy(INPUT.blasthits,upper_margin = 0, lower_margin = 2, remove = c("unwanted_taxon1","unwanted_taxon2"))

# NOTE: The number of OTUs and the upper_limit dictates how many species needs to be classified from kingdom to species, and thus the speed of the function.
# NOTE: The 
# NOTE: In the longer runs a last function to evaluate the classification in relation to the match ID could be included.
#       Of course a best match of 80% means that the species level classification cannot be used. But this is not a part of the script at present.

# Explanation to input
#   INPUT.basthits is the blast-results
#   upper_margin is the margin used for suboptimal hits used for classification - e.g. a margin of 0.5 means that hits of 100% to 99.5% is used og 95% to 94.5%, if the best hit is 100% or 95% respectively.
#   lower_margin: hits down to this margin from the best hit are shown in the output as alternative possibilities, but not used for taxonomic evaluation.
#   remove: a vector of taxa to exclude from the evaluation. Could be e.g. remove = c("uncultured","environmental") to exclude hits with no precise annotation, or names of species known not to be in the study area.
#
# Explanation to output
#   the output is a list with
#   $classified_table : the table with all OTUs classified. One row per OTU
#        this table contains the estimated best classification at all taxonomic levels, based on a weighted score (of the evalue) of the hits in the upper_margin, 
#        each taxonomic level gets a score indicating the agreement on the selected classification at that level..
#        also a string of alternatives and their matches (%) this string includes hits from both upper and lower margin
#   $all_classifications: this is the table used to make the classified_table. It contains all hits above lower_magrin for all OTUs and their classifications (only upper_margin).
#   ...and the input parameters

options(ENTREZ_KEY = "870bca84f5753d51f7cae3f745ff9adb2d09")

# The scripts needs the following packages
library(taxize)
library(dplyr)
library(tidyr)

# EXAMPLE DATA
# In the folder here is a file with some test data, and the resulting output: test_otus.fasta, BLAST_HIT_OUTPUT (result of blast command) and my_classified_otus.txt as the output from the script here.

#This an example to read and format a blast table resulting from the blast command shown above
IDtable=read.csv("BLAST_HIT_OUTPUT",sep='\t',header=F,as.is=TRUE)
#Assign names for columns. Depends on the exact blast command that was excecuted!
names(IDtable) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","qcovs","sgi","sseq","ssciname","staxid")

# Suggestion: Make a small test input to see if it runs.
#IDtable <- IDtable[1:1000,]

# After loading the 4 functions below, this is the command used to classify the data in one go! upper_margin, lower_margin and remove can be adjusted
# We here use an upper margin of 0% to only best hits to evaluate taxonomy, and a lower margin of 2% to include a bit more species in the output to evaluate potential misclassifications.
# Other values can be used. For example an upper_limit of 0.5% to include slightly suboptimal hits in the classification
my_clasified_result <- assign_taxonomy(IDtable,upper_margin = 1, lower_margin = 2, remove = c("uncultured","environmental","N/A"))

#E.g. include evaluation (and taxonomic string) of all hits down to 10% below max
my_clasified_result <- assign_taxonomy(IDtable,upper_margin = 10, lower_margin = 10, remove = c("uncultured","environmental"))

# take a look at the data and save it
my_clasified_result$classified_table
my_clasified_result$all_classifications # Here you can review all scores for all matches per OTU
my_clasified_result$all_classifications_summed # Here you can review all summed scores for all taxa per OTU
write.table(my_clasified_result$classified_table, "my_classified_otus.txt", sep = "\t", quote = F, row.names = F)
write.table(my_clasified_result$adjusted_classified_table, "my_classified_otus_adjusted.txt", sep = "\t", quote = F, row.names = F)

#THESE FOUR FUNCTIONS NEED TO BE LOADED BEFORE RUNNING THE CLASSIFICATION.

#Function1
#Wrapper function using the three main functions - each step can be done manually also...
assign_taxonomy <- function(table,upper_margin=0.5,lower_margin=2, remove = c("uncultured","environmental","N/A"), useDB=TRUE, appendDB=TRUE){
 pf <- prefilter(table, upper_margin, lower_margin, remove)
 gc <- get_classification(pf, useDB, appendDB)
 cf <- evaluate_classification(gc, remove)
 ac <- adjust_classification(cf$taxonon_table)
 result <- list(adjusted_classified_table= ac, classified_table=cf$taxonon_table, all_classifications=cf$all_taxa_table, all_classifications_summed=cf$all_taxa_table_summed,upper=upper_margin, lower=lower_margin, removed=remove)
 return(result)
}

#Function2
#Filter data OTU wise according to upper and lower margin set, and taxa to exclude
prefilter <- function(IDtable, upper_margin=0.5, lower_margin=2, remove = c("")){
 print(paste0("Collapsing taxids and preprocessing")) # make a progressline (indicating the index the loops needs to be
 new_IDtable <- IDtable[0,] # prepare filtered matchlist
 IDtable <- IDtable[!IDtable$staxid == "N/A",]
 ids <- names(table(IDtable$qseqid))
 i=1
 o=length(ids)
 pb <- txtProgressBar(min = 0, max = o, style = 3)
 for (name in ids){
  test <- IDtable[which(IDtable$qseqid == name),] # select all lines for a query
  if (nchar(remove[1])>0){
   test2 <- test
   for (rm in 1:length(remove)){
    test2 <- test2[!grepl(remove[rm], test2$ssciname,ignore.case = TRUE),]
   }
   if (nrow(test2) > 1) {test <- test2}
  }
  max <- max(test$pident)
  upper <- max-upper_margin
  lower <- max-lower_margin
  test <- test[which(test$pident >= lower),] # select all lines for a query
  test$margin <- "lower"
  test[test$pident >= upper,"margin"] <- "upper"
  
  new_IDtable = rbind(new_IDtable,test) # add this row to the filtered IDtable
  i=i+1
  setTxtProgressBar(pb, i)
 }
 close(pb)
 return(new_IDtable)
}

#Function3
# Get full taxonomic path for all hits within the upper limit of each OTU. Identical species are only queried once....
get_classification <- function(IDtable2, useDB=TRUE, appendDB=TRUE, db_path="~/tax_db", tax_db="storedTaxidsRDS"){
 require(taxize)
 all_staxids <- names(table(IDtable2$staxid[IDtable2$margin=="upper"])) # get all taxids for table
 all_classifications <- list() # prepare list for taxize output
 
 if(useDB){
  if(!file_test("-d", db_path)) dir.create(db_path)
  db_file <- file.path(db_path, tax_db)
  if(!file_test("-f", db_file)){
    print("No database file (storedTaxidsRDS) available, classifying all taxids anew")
    useDB=FALSE
   } else {
    storedTaxids <- readRDS(db_file)
    stored_vector <- storedTaxids$taxids %in% all_staxids
    classify_vector <- !all_staxids %in% storedTaxids$taxids
    print(paste0("Using ",sum(stored_vector)," stored classified taxids from database"))
    reused_taxids <- storedTaxids$taxids[stored_vector]
    reused_classifications <- storedTaxids$classifications[stored_vector]
    all_staxids <- all_staxids[classify_vector]
  }
 }
 
 o=length(all_staxids) # number of taxids
 
 print(paste0("Fetching classifications for ",o," taxids from scratch"))
 
 Start_from <- 1 # change if loop needs to be restarted due to time-out
 
 #Get ncbi classification of each entry
 if(o>0){
 pb <- txtProgressBar(min = 0, max = o, style = 3)
 for (cl in Start_from:o){ # the taxize command "classification" can be run on the all_staxids vector in one line, but often there is
  #a timeout command, therefor this loop workaround.
  #restarted from if it quits)
  all_classifications[cl] <- classification(all_staxids[cl], db = "ncbi")
  if(round(cl/10) == cl/10){Sys.sleep(2)}
  setTxtProgressBar(pb, cl)
 }
 close(pb)
 }
 

 if(appendDB){
  db_file <- file.path(db_path, tax_db)
  if(!file_test("-f", db_file)){
   print(paste0("No database file found. Saving all ",length(all_staxids),"classified taxids in a new database (storedTaxidsRDS)"))
   append_tax <- list(taxids=all_staxids,classifications=all_classifications)
   saveRDS(append_tax,db_file)
  } else {
   storedTaxids <- readRDS(db_file)
   append_tax <- list(taxids=c(storedTaxids$taxids,all_staxids),classifications=c(storedTaxids$classifications,all_classifications))
   print(paste0("Appending classification for ",length(all_staxids)," taxids to database"))
   saveRDS(append_tax,db_file)
  }
 }
 
 if(useDB){
  all_staxids <- c(all_staxids,reused_taxids)
  all_classifications <- c(all_classifications,reused_classifications)
 }

 
 #Construct a taxonomic path from each classification
 output <- data.frame(staxid=character(),kingdom=character(), phylum=character(),class=character(),order=character(),family=character(),genus=character(),species=character(), stringsAsFactors=FALSE)
 totalnames <- length(all_staxids)
 pb <- txtProgressBar(min = 0, max = totalnames, style = 3)
 for (curpart in seq(1:totalnames)){
  currenttaxon <- all_classifications[curpart][[1]]
  if (nchar(currenttaxon[1]) > 0) {
   spec <- all_staxids[curpart]
   output[curpart,"kingdom"] <- currenttaxon[which(currenttaxon$rank == "kingdom"),"name"][1]
   output[curpart,"phylum"] <- currenttaxon[which(currenttaxon$rank == "phylum"),"name"][1]
   output[curpart,"class"] <- currenttaxon[which(currenttaxon$rank == "class"),"name"][1]
   output[curpart,"order"] <- currenttaxon[which(currenttaxon$rank == "order"),"name"][1]
   output[curpart,"family"] <- currenttaxon[which(currenttaxon$rank == "family"),"name"][1]
   output[curpart,"genus"] <- currenttaxon[which(currenttaxon$rank == "genus"),"name"][1]
   output[curpart,"species"] <- currenttaxon[which(currenttaxon$rank == "species"),"name"][1]
   output[curpart,"staxid"] <-  spec # add that row to the filtered IDtable
  }
  setTxtProgressBar(pb, curpart)
 }
 close(pb)
 taxonomic_info <- merge(IDtable2,output,by = "staxid", all=TRUE)
 taxonomic_info$species[is.na(taxonomic_info$species)] <- taxonomic_info$ssciname[is.na(taxonomic_info$species)]
 return(taxonomic_info)
}

#Function4
#Function for evaluating the taxonomic assignment of each OTU. All hits within the upper margin are used in the evaluation weithted by thei evalue, so that suboptimal matches has a lower weight. All hits within the lower margin are put into the output (but not used for evaluating classification)
evaluate_classification <- function(classified, remove=c("")){
 require(tidyr)
 require(dplyr)
 ids <- names(table(classified$qseqid))
 i <- 1
 print("Evaluating the classifications")
 pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
 for (name in ids){
  setTxtProgressBar(pb, i)
  test <- classified[which(classified$qseqid == name),]
  test2 <- test %>% filter(margin == "upper")
  if (nchar(remove[1])>0){
   test2x <- test2
   for (rm in 1:length(remove)){
    test2x <- test2x[!grepl(remove[rm], test2x$species,ignore.case = TRUE),]
   }
   if (nrow(test2x) > 1) {test2 <- test2x}
  }
  
  test2$score <- 100*(1/test2$evalue)/sum(1/test2$evalue)  # HER BEREGSES SCOREN FOR ALLE MATCHES PER OTU
  #test2 %>% group_by(species) %>% mutate(n_obs=n())
  
  test4 <- test2 %>% filter(margin == "upper") %>%
   dplyr::select(margin,qseqid,sgi,sseq,staxid,pident,score,qcovs,kingdom,phylum,class,order,family,genus,species) %>% 
   group_by(qseqid,kingdom, phylum,class,order,family,genus,species) %>% 
   mutate(species_score=sum(score)) %>% 
   group_by(qseqid,kingdom, phylum,class,order,family,genus) %>% 
   mutate(genus_score=sum(score)) %>%
   group_by(qseqid,kingdom, phylum,class,order,family) %>% 
   mutate(family_score=sum(score))%>%
   group_by(qseqid,kingdom, phylum,class,order) %>% 
   mutate(order_score=sum(score)) %>%
   group_by(qseqid,kingdom, phylum,class) %>% 
   mutate(class_score=sum(score)) %>%
   group_by(qseqid,kingdom, phylum) %>% 
   mutate(phylum_score=sum(score)) %>%
   group_by(qseqid,kingdom) %>% 
   mutate(kingdom_score=sum(score)) %>% ungroup() %>%
   arrange(-kingdom_score,-phylum_score,-class_score,-order_score,-family_score,-genus_score,-species_score)
  test3 <- test4 %>% slice(1)
  test5 <- test4 %>% distinct(qseqid,sgi,sseq,pident,qcovs,kingdom, phylum,class,order,family,genus,species,kingdom_score,phylum_score,class_score,order_score,family_score,genus_score,species_score) 
  string1 <- test %>% dplyr::group_by(species,pident) %>%  summarize(count=n()) %>% select(species,count,pident) %>% arrange(-pident,-count) %>% t()
  string2 <- toString(unlist(string1))
  test3$alternatives <- string2
  if (i == 1){result <- test3} else{
   result <- rbind(result,test3)
  }
  if (i == 1){result2 <- test2} else{
   result2 <- rbind(result2,test2)
  }
  if (i == 1){result3 <- test5} else{
   result3 <- rbind(result3,test5)
  }
  i=i+1
 }
 close(pb)
 total_result <- list(taxonon_table = result, all_taxa_table=result2, all_taxa_table_summed=result3)
 return(total_result)
}

# Function for adjusting taxonomic annotation. Annotation is adjusted based on the level of match (id_cut). Default settings, applicable for ITS2 data, 
# is 98, 90, 85, 80, 75, 70, 50 for species, genus, family, order, class, phylum, kingdom assignment.
# Assignment is also adjusted for taxonomic agreement among reference database matches. Default threshod scores are 90 for accepting a match at any taxoomic level.

adjust_classification <- function(class_table, id_cut=c(98, 90, 85, 80, 75, 70, 50), uncertainty_levels=c(90,90,90,90,90,90,90)){
 
 my_classifications <- class_table
 
 cutoff_index <- which(my_classifications$pident < id_cut[7])
 reclassify_index <- c("kingdom","phylum","class","order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unidentified"
 my_classifications$species[cutoff_index] <- paste0("unidentified_",my_classifications$phylum[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[6]  & my_classifications$pident >= id_cut[7] & my_classifications$kingdom != "unidentified")
 reclassify_index <- c("phylum","class","order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unidentified"
 my_classifications$species[cutoff_index] <- paste0("unidentified_",my_classifications$phylum[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[5]  & my_classifications$pident >= id_cut[6] & my_classifications$phylum != "unidentified")
 reclassify_index <- c("class","order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unidentified"
 my_classifications$species[cutoff_index] <- paste0("unidentified_",my_classifications$phylum[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[4] & my_classifications$pident >= id_cut[5] & my_classifications$class != "unidentified")
 reclassify_index <- c("order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unidentified"
 my_classifications$species[cutoff_index] <- paste0("unidentified_",my_classifications$class[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[3] & my_classifications$pident >= id_cut[4] & my_classifications$order != "unidentified")
 reclassify_index <- c("family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unidentified"
 my_classifications$species[cutoff_index] <- paste0("unidentified_",my_classifications$order[cutoff_index],"_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[2] & my_classifications$pident >= id_cut[3] & my_classifications$family != "unidentified")
 reclassify_index <- c("genus")
 my_classifications$genus[cutoff_index] <- "unidentified"
 my_classifications$species[cutoff_index] <- paste0("unidentified_",my_classifications$family[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[1] & my_classifications$pident >= id_cut[2] & my_classifications$genus != "unidentified")
 my_classifications$species[cutoff_index] <- paste0("unidentified_",my_classifications$genus[cutoff_index], "_sp")
 
 my_classifications$kingdom[my_classifications$kingdom_score<uncertainty_levels[7]] <- "uncertain"
 my_classifications$phylum[my_classifications$phylum_score<uncertainty_levels[6]] <- "uncertain"
 my_classifications$class[my_classifications$class_score<uncertainty_levels[5]] <- "uncertain"
 my_classifications$order[my_classifications$order_score<uncertainty_levels[4]] <- "uncertain"
 my_classifications$family[my_classifications$family_score<uncertainty_levels[3]] <- "uncertain"
 my_classifications$genus[my_classifications$genus_score<uncertainty_levels[2]] <- "uncertain"
 
 print("Adjusting the classifications")
 pb <- txtProgressBar(min = 0, max = nrow(my_classifications), style = 3)
 
 for(i in 1:nrow(my_classifications)){
  setTxtProgressBar(pb, i)
  if (my_classifications[i,"species_score"] < uncertainty_levels[1]){
   if (!my_classifications[i,"genus"] %in% c("uncertain", "unidentified", NA)){
    my_classifications[i,"species"] <- paste0("uncertain_",my_classifications[i,"genus"],"_sp")
   } else if (!my_classifications[i,"family"] %in% c("uncertain", "unidentified", NA)){
    my_classifications[i,"species"] <- paste0("uncertain_",my_classifications[i,"family"],"_sp")
   } else if (!my_classifications[i,"order"] %in% c("uncertain", "unidentified", NA)){
    my_classifications[i,"species"] <- paste0("uncertain_",my_classifications[i,"order"],"_sp")
   } else if (!my_classifications[i,"class"] %in% c("uncertain", "unidentified", NA)){
    my_classifications[i,"species"] <- paste0("uncertain_",my_classifications[i,"class"],"_sp")
   } else if (!my_classifications[i,"phylum"] %in% c("uncertain", "unidentified", NA)){
    my_classifications[i,"species"] <- paste0("uncertain_",my_classifications[i,"phylum"],"_sp")
   } else if (!my_classifications[i,"kingdom"] %in% c("uncertain", "unidentified", NA)){
    my_classifications[i,"species"] <- paste0("uncertain_",my_classifications[i,"kingdom"],"_sp")
   }
  }
 }
 close(pb)
 return(my_classifications)
}

