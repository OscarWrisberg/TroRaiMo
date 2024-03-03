############################
# Taxon matching taxonomy to WCVP #
############################
#Setting the dirname

###################################################################################
#################### Settings for local testing of the script #####################
###################################################################################

###setting the wd
#setwd("/home/au543206/GenomeDK/Trf_models/workflow/01_distribution_data/04_common_format")
# input_file <- "gbif_common_format.rds"
# wcvp_input_file <- "wcvp_names_apg_aligned.rds"
# output_file <- "gbif_taxon_matched.rds"

###################################################################################
###################################################################################

# Loading packages
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table","dplyr","tidyr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
#####################################################################################

DB.name="GBIF"

args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1]) # here you define the name of your input file
wcvp_input_file <- as.character(args[2]) # Here you define the name of the WCVP file you would like to use.
output_file <- as.character(args[3]) # Here you define the name of the output file

###############################################################################################################################
################################################### READ + ADJUST WCVP DATA ###################################################
###############################################################################################################################

# Load data
cat("Started to read Kew inputfile at: ", format(Sys.time(),"%H:%M:%S"), "\n")
wc_all <- readRDS(wcvp_input_file)
cat("Finished reading Kew inputfile at: ", format(Sys.time(),"%H:%M:%S"), "\n")

# add new column "tax_comb" to unify taxon ranks and infraspecific ranks column
wc_all$tax_comb <- wc_all$infraspecific_rank # transfer content
wc_all$tax_comb <- gsub("\\.|,", "", wc_all$tax_comb) # clean up, remove dots and commas
empty_ranks <- which(wc_all$tax_comb == "" | is.na(wc_all$tax_comb)) # get empty infraspec ranks
wc_all$tax_comb[empty_ranks] <- as.character(wc_all$taxon_rank[empty_ranks]) # Here it adds the information from the taxon_rank column

# remove capitals
wc_all$tax_comb <- gsub("Species", "species", wc_all$tax_comb)
wc_all$tax_comb <- gsub("Genus", "genus", wc_all$tax_comb)

# # replace empty cells with NA
wc_all[wc_all == ""] <- NA


# fix hybrid notation
#wc_all$genus_hybrid <- as.character(wc_all$genus_hybrid)
#wc_all$species_hybrid <- as.character(wc_all$species_hybrid)
wc_all$genus_hybrid[!is.na(wc_all$genus_hybrid)] <- "x"
wc_all$species_hybrid[!is.na(wc_all$species_hybrid)] <- "x"

# subset the dataset with key columns and remove all NA accepted IDs
wc_all_sub <- wc_all %>%
  dplyr::select(tax_comb, taxon_status, family.apg, genus, genus_hybrid, species, 
                species_hybrid, infraspecies, taxon_authors, accepted_plant_name_id) %>% 
  unique() %>%
  filter(!is.na(accepted_plant_name_id))

# rename columns
names(wc_all_sub)[names(wc_all_sub) %in% c("tax_comb", "infraspecies", "taxon_authors")] <-
  c("taxon_rank", "infra_name", "author")



##################################################################################
##############################  GBIF  ############################################
##################################################################################
if (DB.name == "GBIF") {

  cat("Started to read GBIF inputfile at: ", format(Sys.time(), "%H:%M:%S"), "\n")
  input <- readRDS(input_file)
  cat("Finished reading GBIF inputfile at: ", format(Sys.time(), "%H:%M:%S"), "\n")

  # RENAMING family column
  input$family.apg <- input$family
  input <- input %>%
    dplyr::select(-family)
  
  # rename
  input$taxon_rank <- gsub("variety", "var", input$taxon_rank, fixed = TRUE)
  input$taxon_rank[which(is.na(input$taxon_rank))] <- "undefined"
  
  # set data classes right: all factors need to be changed to characters
  classes <- as.vector(unlist(lapply(input, class)))
  input[, which(classes == "factor")] <- lapply(input[, which(classes == "factor")], as.character)
  
  # exclude taxa without extracted genus
  if(any(is.na(input$genus))){
    input <- input[-which(is.na(input$genus)),]      
  }
  
  dataset <- input
  rm(input)
}  

##########################################################################################################################################################################


#### TAXONOMY MATCHING - Shared Functionality #################################

###############################################################################
######################## STRICT MATCHING, no author ###########################
###############################################################################
# In this section of the code we are matching the taxon names from the GBIF data frame to the WCVP data frame.
# The strict matching is done by matching the taxon_rank, family.apg, genus, genus_hybrid, species, species_hybrid and infra_name columns.
# If all of these columns match then the taxon name is a strict match.
# The strict matches are then looked through to see if there are any multimatches.
# These are assigned to their own data frame.
# at the end we check if all the species in the dataset are accounted for in either the multimatches, no matches or strict matches

# Merging occurrences and WCVP
# This is SUPER important as this is deemed a strict match if the fields match
# Why does it result in all the accepted_plant_name_id's being NA, apparantly because there are no matching columns...
res <- merge(dataset, wc_all_sub, all.x = TRUE,
             by=c("taxon_rank", 
                 "genus",
                 "family.apg",
                 "genus_hybrid",
                 "species",
                 "species_hybrid",
                 "infra_name"
              ))


# get rid of duplicate columns
res.tmp <- res %>%
  dplyr::select(-author.x, -author.y) %>%
  unique() %>%
  mutate(match_type=ifelse(is.na(accepted_plant_name_id), NA, "strict match"))


# get IDs that have no match
no_match <- res.tmp %>%
  filter(is.na(accepted_plant_name_id))
match <- res.tmp %>% 
  filter(!is.na(accepted_plant_name_id))

# Report on how big a proportion of the total data can be strictly matched
cat(" \n Proportion of strict matches is ", length(match$id)/length(dataset$id), "\n \n")

# get IDs that have multiple matches in the strict matching
mm_strict_id <- unique(match$id[duplicated(match$id)])

# Report on how big a proportion of the strict matches are multimatches
cat(paste0("\n \n Proportion of strict matches that are multimatches is ", length(mm_strict_id)/length(match$id), "\n \n"))

# get single matches
done <- match %>%   
  filter(!(id %in% mm_strict_id)) %>% # filtering out the multimatches
  dplyr::select(id, accepted_plant_name_id, match_type) # selecting the columns we want

# Report on the proportion of single matches
cat(paste0("\n \n Proportion of strict matches that are single matches is ", length(done$id)/length(match$id), "\n \n"))

# get multimatches
mm <- res %>% 
  filter(id %in% mm_strict_id) %>% 
  arrange(id)

# Remove duplicates from multimatches
mm <- distinct(mm)

# remove taxa from no_match if they occur in done or mm
no_match <- no_match %>% 
  filter(!(no_match$id %in% done$id))
no_match <- no_match %>% 
  filter(!(no_match$id %in% mm$id))

# last overall check
try(if(all(unique(c(mm$id, no_match$id, done$id)) %in% dataset$id)){
  print("all IDs represented!")
}else{
  stop("NOT all IDs represented!")
})




##############################################################################################
#################################### MULTIMATCHES ############################################
##############################################################################################
# In this piece of code we check the multimatches to see if they can be resolved.
# We construct a multi_match_checker function which looks through the multimatched ID's
  # This multi_match_checker checks if all the matches for the ID, point to the same accepted ID.
  # If they do then the ID is a resolved multimatch
  # otherwise it checks if there is only one of the matches which also has a matching author.
  # If there is then the ID is a resolved multimatch
  # If there is no author or more than one match with missing author, then the ID is an unresolved multimatch

# At the end, we update the done data frame with the resolved multimatches
# and the unresolved multimatches are saved to their own dataframe.

# simplify author names: remove all spaces, dots and numbers
mm$author.x <- gsub(" |\\.|[0-9]*|,", "", mm$author.x) # author names from dataset
mm$author.y <- gsub(" |\\.|[0-9]*|,", "", mm$author.y) # author names from WCVP

# assign matching authors status
mm$match_type <- ifelse(mm$author.x == mm$author.y, "matching authors", "no matching authors")

#### define a function to check multiple matches ####
# This function checks how 
multi_match_checker <- function(mm_ids, subdata){
  valid_matches <- c("matching authors", "strict match no family", "no infra (same same)")
  one_match <- NULL
  more_match <- NULL
  zero_match <- NULL
  resolved_mm <- NULL
  unresolved_mm <- NULL
  total_mm_ids <- length(sort(unique(mm_ids)))
  sorted_unique_ids <- sort(unique(mm_ids))

  cat("Total number of unique MM ID's supplied is ",total_mm_ids, "\n ")
  
  # loop over all IDs
  for (i in seq_along(sorted_unique_ids)){
    temp <- subdata[subdata$id==sorted_unique_ids[i],]
     
    # print progress
    if(!i%%1000) cat("Percentage done",format(round((i/total_mm_ids)*100,3), nsmall = 3)," at ",format(Sys.time(),'%H:%M:%S'), "\n")

    # all entries point to the same accepted ID
    if(length(unique(temp$accepted_plant_name_id)) == 1){
      one_match <- c(one_match, sorted_unique_ids[i])
      resolved_mm <- rbind(resolved_mm, temp)
      # entries point to different accepted IDs
    }else{
      # one entry with matching authors
      if(sum(temp$match_type %in% valid_matches, na.rm = TRUE) == 1){
        one_match <- c(one_match, sorted_unique_ids[i])
        resolved_mm <- rbind(resolved_mm, temp[which(temp$match_type %in% valid_matches),])
      }
      # > 1 entries with matching authors
      if(sum(temp$match_type %in% valid_matches, na.rm=TRUE) > 1){
        more_match <- c(more_match, sorted_unique_ids[i])
        unresolved_mm <- rbind(unresolved_mm, temp)
      }
      # no matches because missing authors
      if(all(is.na(temp$match_type))){
        zero_match <- c(zero_match, sorted_unique_ids[i])
        unresolved_mm <- rbind(unresolved_mm, temp)
      }
    }
  }
  results <- list(
    "one_match" = one_match,
    "more_match" = more_match,
    "zero_match" = zero_match,
    "resolved_mm" = resolved_mm,
    "unresolved_mm" = unresolved_mm
  )
  # Print in terminal the proportion of zero matches, one matches and more matches.
  cat("Proportion unmatchable due to missing authors ", length(zero_match)/total_mm_ids, "\n")
  cat("Proportion of matches with one accepted ID ", length(one_match)/total_mm_ids, "\n")
  cat("Proportion of matches where there are more than one valid match ", length(more_match)/total_mm_ids, "\n")
  #cat("Proportion of resolved matches is ", length(resolved_mm$id)/total_mm_ids, "\n")
  #cat("Proportion of unresolved matches is ", length(unresolved_mm$id)/total_mm_ids, "\n")
  return(results)
}

# Writing out the dimensions of mm
#cat(paste0("Length of mm_strict_id is: ", length(mm_strict_id), " \n"))
#cat(paste0(" Dimensions of mm is: ", dim(mm), " \n"))

cat(paste0("Starting multi_match_checker at 1st time at: ", format(Sys.time(),'%H:%M:%S'), "\n"))

#This function takes around 40 hours to run on a cluster with 75 gigs of ram and 15 cores.
# Here mm_strict_id is the unique id's from the mm data frame where the id's are duplicated
# and mm is the dataframe of the multimatches
mm_results <- multi_match_checker(mm_strict_id, mm)

cat(paste0("Ended multi_match_checker at: \n", Sys.time(), "\n"))

# if unresolved multimatces are present save them into a separate file
if(!is.null(mm_results$unresolved_mm)){ 
wcp_conflicts1 <- mm_results$unresolved_mm %>% 
  dplyr::filter(id %in% mm_results$more_match)
}

# update done data set with resolved entries
if(length(mm_results$resolved_mm > 0)){
done <- rbind(done, mm_results$resolved_mm %>% 
                dplyr::select(id, accepted_plant_name_id, match_type))

}






##############################################################################
######################### NO STRICT MATCHES ##################################
##############################################################################
# Here we are checking for matches for hybrid species.
# Firstly we take the res dataframe and select only the id's which are in the no_match dataframe.
# We then create a subset of the no_match dataframe by filtering out entries where the infraspecific name is the same as the species name and if the id is not in the done dataframe.

#This subset is then merged with the wc_all_sub dataframe by the columns family.apg, genus, genus_hybrid, species, species_hybrid and author.
# If the accepted_plant_name_id is NA then the match_type is NA otherwise it is "no infra (same same)"
# This merged dataframe is then saved and we will then check for multimatches in this dataframe.

# In order to check for multimatches we first find all the id's which are duplicated in the res3 dataframe.
# We then create a subset of the res3 dataframe by filtering out the id's which are in the mm_strict_id dataframe and where the match_type is not NA.
# This subset is then supplied to the multi_match_checker function from above and the resolved multimatches are then saved in the done dataframe.

## no infraspecific match: if infra == species name → match on species level

no_match <- res %>% # Here we are taking the res data frame from above
  dplyr::filter(res$id %in% no_match$id) %>% # Here we are filtering out the no_match id's from the res data frame
  dplyr::select(-author.y, -taxon_status, -accepted_plant_name_id) %>% # Here we are removing the author.y, taxon_status and accepted_plant_name_id columns
  unique() %>% #
  dplyr::rename(author = author.x) # Here we are renaming the author.x column to author

## exclude taxa that have been matched to avoid double assignment
nms <- no_match %>% 
  filter(infra_name == species & !(id %in% done$id))

# merge with WCVP if match_type is NA otherwise write that match_type is "no infra (same same)"
res3 <- merge(nms, wc_all_sub, all.x=TRUE, # 
              by=c("family.apg", "genus", "genus_hybrid", 
                   "species", "species_hybrid", "author")) %>% 
  mutate(match_type=ifelse(is.na(accepted_plant_name_id), NA, "no infra (same same)"))

## get species ID with multiple matches
mm_strict_id <- unique(res3$id[duplicated(res3$id)]) # Here we are getting the unique id's from the res3 data frame where the id's are duplicated
res3_sub <- res3 %>% 
# Here we are filtering out the id's from the res3 data frame where the id's are in the mm_strict_id data frame and where the match_type is not NA
  filter(!(id %in% mm_strict_id) & !is.na(match_type))

# attach single matches to the resolved dataframe
done <- rbind(done, res3_sub %>% 
                dplyr::select(id, accepted_plant_name_id, match_type))

# generate multimatches dataset
res3_mm <- res3 %>% 
  filter(id %in% mm_strict_id)

mm_results3 <- multi_match_checker(mm_strict_id, res3_mm)

# attach the resolved multimaches to done
if(!all(lapply(mm_results3, length)==0)){
  done <- rbind(done, mm_results3$resolved_mm %>% 
                  dplyr::select(id, accepted_plant_name_id, match_type))
}






#####################################################################################
######################### no family, but authors ####################################
#####################################################################################
# Here we are checking for matches for hybrid species.
# but this time were are checking for matches where the family is not the same but the author is the same.
# Firstly we take the res dataframe and select only the id's which are in the no_match dataframe.


## exclude taxa that have been matched already to avoid double assignment
no_match <- no_match %>% 
  filter(!(id %in% done$id))

res2 <- merge(no_match, wc_all_sub, all.x=TRUE, 
              by=c("taxon_rank", "genus", "genus_hybrid",
                   "species", "species_hybrid", "infra_name", "author")) %>% 
  mutate(match_type=ifelse(is.na(accepted_plant_name_id), NA, "strict match no family"))

## get species ID with multiple matches
mm_strict_id <- unique(res2$id[duplicated(res2$id)])

## attach single matches to the resolved dataframe
res2_sub <- res2 %>% 
  filter(!(res2$id %in% mm_strict_id) & !is.na(match_type))
done <- rbind(done, res2_sub %>% 
                dplyr::select(id, accepted_plant_name_id, match_type))

# generate multimatches dataset, altough differtent families
res2_mm <- res2 %>% 
  filter(id %in% mm_strict_id)

if(length(mm_strict_id)>0){
  mm_results2 <- multi_match_checker(mm_strict_id, res2_mm)
  
  # attach the resolved multimaches to done
  done <- rbind(done, mm_results2$resolved_mm %>% 
                  dplyr::select(id, accepted_plant_name_id, match_type))
}else{
  mm_results2 <- NULL
}

####################################################################################
############## final check and accepted ID assignment ##############################
####################################################################################
# This part checks if there are any ID's in the done dataframe which are duplicated.
# it then removes entries which have duplicated ID's in the done dataframe and where the match type is "no matching authors"


# no double assignments? should all be 1
table(tapply(done$accepted_plant_name_id, done$id, function(x)length(unique(x))), useNA="ifany")


# remove NAs + no matching authors from converging multimatch merge
dupli_match_types <- done$id[which(duplicated(done$id))]
done <- done %>% 
  filter(!is.na(match_type) | !is.na(id)) %>%
  unique()
done <- done[-which(done$id %in% dupli_match_types & done$match_type=="no matching authors"),]

# check for double assignments
if(any(duplicated(done$id))){
  dupl.id <- done$id[duplicated(done$id)]
  print("double assignment happend - check double_assignments for cases")
  double_assignments <- done %>% filter(id %in% dupl.id)
}else print("No double assignments. All is well!")

# attach done data frame to input data set
fin <- merge(dataset, done, by="id", all.x=TRUE)

# combine unresolved into one file
unsolved <- bind_rows(mm_results$unresolved_mm, mm_results2$unresolved_mm, mm_results3$unresolved_mm)

# assign multimatch 
if(any(fin$id %in% unsolved$id)){
  fin$match_type[fin$id %in% unsolved$id] <- "multimatch"
}

wcsp.acc.name <- wc_all %>% 
  dplyr::select(taxon_name, taxon_status, accepted_plant_name_id) %>% 
  filter(taxon_status=="Accepted")

wcsp.w.acc.name <- left_join(fin, wcsp.acc.name, by="accepted_plant_name_id")


#######################################################################################
########################### SPECIES LEVEL MATCHING ####################################
#######################################################################################

# matching logic: get accepted ID for all subspecies entries 
# --> get genus + species name from those accepted IDs names
# --> create another dataframe from wcsp only containing species level entries
# --> match the subspecies-level dataframe with the species-level dataframe based on genus and species name
# --> attach the species-level accepted ID to the final dataframe

## subset wcsp to species present in input data
ids <- unique(fin$accepted_plant_name_id)
wc_sub <- wc_all[wc_all$plant_name_id %in% ids,]

## subset to taxon rank < species
wc_subspecies <- wc_sub[!wc_sub$tax_comb %in% c("species", "genus"),]

# get species-level entries from wcsp. (do they have to be accepted?)
wc_species <- wc_all[wc_all$tax_comb %in% c("species", "genus") & wc_all$taxon_status=="Accepted",]

# match subsepcies with species based on which have the same genus + species entry
## plant name id from subspecies will be used to match to the fin dataframe using accepted ID
species_match <- merge(wc_subspecies[,c("plant_name_id", "species", "genus")],
                       wc_species[,c("accepted_plant_name_id", "species", "genus")],
                       by=c("species", "genus"),
                       all.x=TRUE)

# there are XXX multimatches
paste0("There are ",(nrow(species_match)-nrow(wc_subspecies)), " multimatches")

species_match$elevated_to_species_id <- species_match$accepted_plant_name_id

## get species ID with multiple matches
multimatch_id <- unique(species_match$plant_name_id[which(duplicated(species_match$plant_name_id))])
species_mms <- species_match[species_match$plant_name_id %in% multimatch_id,]

species_match_sub <- species_match[!species_match$plant_name_id %in% multimatch_id,]
species_match_sub <- species_match_sub[!is.na(species_match_sub$accepted_plant_name_id),]

#add to fin
fin_species_match <- merge(fin, species_match_sub[,c("plant_name_id", "elevated_to_species_id")],
                           by.x="accepted_plant_name_id",
                           by.y="plant_name_id", all.x=TRUE)

# save output
saveRDS(fin_species_match, file=output_file)


# Note that some species will not have species level ID assigned, that is
# because the accepted plant name ID is at species level already (example:
# "Equisetum variegatum var. alaskanum") If the original taxon name input was
# species level, and it still gets an elevated to species level ID assigned,
# that means the species name points to an accepted subspecies, which again was
# matched with its parent species name  (see for example "Centaurea asperula")

