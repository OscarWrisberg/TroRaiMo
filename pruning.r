# Load the required packages
library(ape)
library(data.table)
library(dplyr)
library(phytools)
library(stringdist)

###########################################################################################################################
# Settings for running this script locally.
setwd("/home/au543206/GenomeDK/Trf_models/data") # Set the working directory when local
# setwd("/home/owrisberg/Trf_models/data") # Set working directory when remove
wcvp <- readRDS("../workflow/02_adding_orders/wcvp_names_apg_aligned.rds")  # Read the WCVP names file into a data frame
tree <- read.tree("GBMB.tre") # Read the GBMB tree
output_file_tree <- "GBMB_pruned.tre" # Set the name of the output file
###########################################################################################################################

########################################
### Running the script on the cluster ##
########################################

# # Command line arguments
# input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
# input_file_wcvp <- commandArgs(trailingOnly = TRUE)[2]
# output_file_tree <- commandArgs(trailingOnly = TRUE)[3]

# # Read the WCVP names file into a data frame
# cat("Opening ", input_file_wcvp, "\n")  
# wcvp <- readRDS(input_file_wcvp)

# # Read the GBMB tree
# cat("Opening ", input_file_tree, "\n")
# tree <- read.tree(input_file_tree)

tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

# Hardcoding in to drop the stupid Passiflora hybrid cultivar which fuzzy matches to Passiflora hybrida.
tree <- drop.tip(tree, "Passiflora hybrid cultivar")

# Get tip names from the tree
tip_names <- tree$tip.label

# Find matching and non-matching tips
matching_tips <- tip_names[tip_names %in% wcvp$taxon_name] # 76935 tips are matching
not_matching_tips <- tip_names[!(tip_names %in% wcvp$taxon_name)] # Only 2939 tips are not matching
matching_tips_accepted <- tip_names[tip_names %in% wcvp$taxon_name[wcvp$taxon_status == "Accepted"]] # 66166 out of the matching tips are accepted
matching_tips_not_accepted <- tip_names[tip_names %in% wcvp$taxon_name[wcvp$taxon_status != "Accepted"]] # 14355 out of the matching tips are not accepted

# Finding length of matching tips
cat("Length of matching tips ", length(matching_tips), "\n")
cat("Length of non-matching tips ", length(not_matching_tips), "\n")
cat("Length of matching tips accepted ", length(matching_tips_accepted), "\n")
cat("Length of matching tips not accepted ", length(matching_tips_not_accepted), "\n")
# Create a data frame with tip names and families

# Looking at the not_matching_tips I need to loop through them in order to check if they are indeed unmatchable or if they can be matched.
# I could use the taxonomy_matcher to find the correct name for all the tips in the tree, but I think this would be too time consuming.
# Instead I am just going to loop through the names and see if I can find a match in the WCVP file using the grep function.

# # Initialize empty vectorsq
# not_matchable_tips <- character(0)
# matchable_tips <- character(0)
# match_name <- character(0)

# # Loop through the not_matching_tips
# for (i in seq_along(not_matching_tips)) {
#   # Writing a progress bar
#   if (!i %% 50) cat("Percentage done", format(round((i / length(not_matching_tips)) * 100, 2), nsmall = 2), " at ", format(Sys.time(), '%H:%M:%S'), "\n")

#   # Check for exact match
#   if (not_matching_tips[i] %in% wcvp$taxon_name) {
#     #cat("Exact match found for ", not_matching_tips[i], "\n")
#     matchable_tips <- c(matchable_tips, not_matching_tips[i])
#     match_name <- c(match_name, not_matching_tips[i])
#   } else {
#     # Check for matches with one substitution, insertion, or deletion
#     potential_matches <- stringdist::stringdistmatrix(not_matching_tips[i], wcvp$taxon_name, method = "lv", useNames = TRUE)
    
#     # Filter for matches with distance 1
#     matches_with_distance_one <- wcvp$taxon_name[potential_matches <= 1]
    
#     # Check if any matches were found
#     if (length(matches_with_distance_one) > 0 & length(matches_with_distance_one) ==) {
#       #cat("Match found for ", not_matching_tips[i], " with one substitution, insertion, or deletion\n")
# 	  if (length(matches_with_distance_one) == 1) {
#       #cat("Match is ", matches_with_distance_one[1], "\n")
#       matchable_tips <- c(matchable_tips, not_matching_tips[i])
#       match_name <- c(match_name, matches_with_distance_one[1])
# 	  } else {
# 		  cat("Multiple matches found for ", not_matching_tips[i], "\n")
# 		  cat("Matches are ", matches_with_distance_one, "\n")
# 		  not_matchable_tips <- c(not_matchable_tips, not_matching_tips[i])
#     } 
# 	} else {
#       # No match found
#       #cat("No match found for ", not_matching_tips[i], "\n")
#       not_matchable_tips <- c(not_matchable_tips, not_matching_tips[i])
#     }
#   }
# }

# # Saving vectors as RDS files
# saveRDS(not_matchable_tips, "not_matchable_tips_1.rds") 
# saveRDS(matchable_tips, "matchable_tips_1.rds") 
# saveRDS(match_name, "match_name_1.rds")

# Loading RDS files
not_matchable_tips <- readRDS("not_matchable_tips_1.rds")
matchable_tips <- readRDS("matchable_tips_1.rds")
match_name <- readRDS("match_name_1.rds")


cat("Are all the matchable tips found in the tree? ", all(matchable_tips %in% tree$tip.label), " and are all the match_names found in the wcvp$taxon_name", all(match_name %in% wcvp$taxon_name), "\n")
# Renaming tips in the tree with the matched tip names
tree$tip.label[which(tree$tip.label %in% matchable_tips)] <- match_name # This works because both of them is ordered

# Are all species accounted for?
cat("Have all the match_name's been incorporated into the tip.labels: ",all(match_name %in% tree$tip.label))
cat("Are all species accounted for: ", all(not_matching_tips %in% not_matchable_tips | not_matching_tips %in% matchable_tips), "\n")

###############################################################################################################################################
# Can I loop through the unmatched tips and find the tips where we have a subsp or a variety? and where we can find the parent species.
# looping through the not_matchable_tips which when split by " " gives me a vector longer than 2.
###############################################################################################################################################

# # # Initialize empty vectors
# split_not_matchable_tips <- character(0)
# split_matchable_tips <- character(0)
# split_match_name <- character(0)
# split_multi_match <- character(0)

# # Loop through the not_matchable_tips
# for (i in seq_along(not_matchable_tips)) {
#   # Writing a progress bar
#   if (!i %% 50) cat("Percentage done", format(round((i / length(not_matchable_tips)) * 100, 2), nsmall = 2), " at ", format(Sys.time(), '%H:%M:%S'), "\n")

#   # Split the tip by space
#   tip_elements <- strsplit(not_matchable_tips[i], " ")[[1]]
  
#   # If the split list is only of length 2, then we can go to next tip
#   if(length(tip_elements) <= 2){
#     split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
# 	next
#   }

#   if( any(tip_elements == "x" | tip_elements == "X") ){
#     # If the tip contains an x or an X then we can go to next tip
#     #cat("Cannot match ", not_matchable_tips[i], " because it contains an x or an X \n")
#     split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
#     next
#   } 

#   if ( any(tip_elements == "Sp." | tip_elements == "sp.")){
#     # If the tip contains an Sp. or an sp. then we can go to next tip
#     #cat("Cannot match ", not_matchable_tips[i], " because it contains an Sp. or an sp. \n")
#     split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
#     next
#   }

#   # Otherwise we can try to find a match for the first two elements
#   if (length(tip_elements) > 2) {
#     # Take the first two elements
#     tip_to_search <- paste(tip_elements[1:2], collapse = " ")
    
#     # Check for exact match
#     if (tip_to_search %in% wcvp$taxon_name) {
#       #cat("Exact match found for ", not_matchable_tips[i], "\n")
#       split_matchable_tips <- c(split_matchable_tips, not_matchable_tips[i])
#       split_match_name <- c(split_match_name, tip_to_search)
#     } else {
#       # Check for matches with one substitution, insertion, or deletion
#       potential_matches <- stringdist::stringdistmatrix(tip_to_search, wcvp$taxon_name, method = "lv", useNames = TRUE)
      
#       # Filter for matches with distance 1
#       matches_with_distance_one <- wcvp$taxon_name[potential_matches <= 1]
      
#       # Check if any matches were found
#       if (length(matches_with_distance_one) > 0 & length(matches_with_distance_one) == 1) {
#         #cat("Match found for ", not_matchable_tips[i], " with one substitution, insertion, or deletion\n")
#         #cat("Match is ", matches_with_distance_one[1], "\n")
#         split_matchable_tips <- c(split_matchable_tips, not_matchable_tips[i])
#         split_match_name <- c(split_match_name, matches_with_distance_one[1])
#       }

#        if(length(matches_with_distance_one) > 0 & length(matches_with_distance_one) > 1) {
#         #cat("Multiple matches found for ", not_matchable_tips[i], "\n")
#         #cat("Matches are ", matches_with_distance_one, "\n")
#         split_multi_match <- c(split_multi_match, not_matchable_tips[i])

#       } else {
#         # No match found
#         #cat("No match found for ", not_matchable_tips[i], "\n")
#         split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
#       }
#     }
#   }
# }

# #Saving vectors as RDS files for easy loading if I need to rerun the script
# cat("Saving Rds files \n")
# saveRDS(split_not_matchable_tips, "split_not_matchable_tips.rds")
# saveRDS(split_matchable_tips, "split_matchable_tips.rds")
# saveRDS(split_match_name, "split_match_name.rds")
# saveRDS(split_multi_match, "split_multi_match.rds")
# cat("Rds files saved \n")

## Loading RDS files
split_not_matchable_tips <- readRDS("split_not_matchable_tips.rds")
split_matchable_tips <- readRDS("split_matchable_tips.rds")
split_match_name <- readRDS("split_match_name.rds")
split_multi_match <- readRDS("split_multi_match.rds")


# Are all species accounted for?
cat("Are all species accounted for: ", all(not_matchable_tips %in% split_not_matchable_tips | not_matchable_tips %in% split_matchable_tips| not_matchable_tips %in% split_multi_match), "\n")

# And now we can again rename the tips based on the matches we found
cat("Number of tips in split_,matchable_tips is: ", length(split_matchable_tips), "\n")
cat("Number of tips in split_match_name is: ", length(split_match_name), "\n")
cat("Number of tips in tree which are found in matchable names are: ", length(tree$tip.label[which(tree$tip.label %in% split_matchable_tips)]))

# # Making a dataframe with the split_matchable_tips and the split_match_name
# split_match_tips_names <- as.data.frame(cbind(split_matchable_tips, split_match_name))

# # Removing tips which are not in the tree.
# for (i in seq_along(split_match_tips_names[which(!split_match_tips_names$split_matchable_tips %in% tree$tip.label),1])){
#   cat("Removing ",split_match_tips_names[which(!split_match_tips_names$split_matchable_tips %in% tree$tip.label),1][i], " from the matchable tips because it is not in the tree \n")
#   split_match_tips_names <- split_match_tips_names[-which(split_match_tips_names$split_matchable_tips == split_match_tips_names[which(!split_match_tips_names$split_matchable_tips %in% tree$tip.label),1][i]),]
# }

# # Renaming tips 
# tree$tip.label[which(tree$tip.label %in% split_match_tips_names$split_matchable_tips)] <- split_match_tips_names$split_match_name 

# Removing the row from split_matchable_tips and split_match_name if the tip is not in the tree
split_match_name <- split_match_name[which(split_matchable_tips %in% tree$tip.label)]
split_matchable_tips <- split_matchable_tips[which(split_matchable_tips %in% tree$tip.label)]

# Renaming tips in the tree with the matched tip names
tree$tip.label[which(tree$tip.label %in% split_matchable_tips)] <- split_match_name # This works because both of them is ordered


# Are all species accounted for?
cat("Have all the split_match_name's been incorporated into the tip.labels: ",all(split_match_name %in% tree$tip.label))

# Remove tips that are still not matched
cat("Removing", length(split_not_matchable_tips) ,"tips that are still not matched \n")
tree <- drop.tip(tree, split_not_matchable_tips)

cat("Removing", length(split_multi_match), "split matches with multiple matches \n")
tree <- drop.tip(tree, split_multi_match)

# Are all the tips in the tree now found in the WCVP file?
cat("Are all the tips in the tree now found in the WCVP file? ", all(tree$tip.label %in% wcvp$taxon_name), "\n")

cat("We managed to find matches for ", length(matchable_tips) + length(split_matchable_tips), " tips out of ", length(not_matching_tips), "not matching tips \n")
cat("but we had to drop ", length(split_not_matchable_tips) + length(split_multi_match), " tips \n")

# Updating the number of species in the tree
# Get tip names from the tree
tip_names <- tree$tip.label
# Find matching and non-matching tips
matching_tips <- tip_names[tip_names %in% wcvp$taxon_name] # 76935 tips are matching
not_matching_tips <- tip_names[!(tip_names %in% wcvp$taxon_name)] # Only 2939 tips are not matching
matching_tips_accepted <- tip_names[tip_names %in% wcvp$taxon_name[wcvp$taxon_status == "Accepted"]] # 66166 out of the matching tips are accepted
matching_tips_not_accepted <- tip_names[tip_names %in% wcvp$taxon_name[wcvp$taxon_status != "Accepted"]] # 14355 out of the matching tips are not accepted

# Finding length of matching tips
cat("Length of matching tips ", length(matching_tips), "\n")
cat("Length of non-matching tips ", length(not_matching_tips), "\n")
cat("Length of matching tips accepted ", length(matching_tips_accepted), "\n")
cat("Length of matching tips not accepted ", length(matching_tips_not_accepted), "\n")

#End of dealing with tips which were not matching to the WCVP file.
# Now all tips in the tree should be found in the WCVP file.
# This means we now need to start dealing with the species in the tree which are not Accepted in the WCVP file.
# And we need to deal with the tips in the tree which are not species but something else (subsp, var, etc.)


############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
# Now we start dealing with the tips in the tree which are found in the wcvp but are not accepted.
# One of the challenges with using the wcvp and taxon_names is that taxon_names are not unique. and therefore they can be both Accepted and Synonyms at the same time.
# In order to find out how many of the species in my tree which are found as both an accepted and a synonym I should be able to use a subset of the wcvp which only
# Contains the accepted species and then count how many of my tips does not appear in this subset


# Now I need to create a version of the wcvp which have only accepted names.
wcvp_accepted <- wcvp[which(wcvp$taxon_status == "Accepted"),]

new_matching_tips_not_accepted <- tree$tip.label[which(!tree$tip.label %in% wcvp_accepted$taxon_name)] # they are all in the old matching tips not accepted.
cat("There are ",length(matching_tips_not_accepted)-length(new_matching_tips_not_accepted) , " tips in the tree which are both Accepted and Synonym due to different authors. \n")
cat("This means there are ",length(new_matching_tips_not_accepted)," tips in the tree which are not accepted at all \n") # still 11044 tips are not accepted, which is down from 14355. This means 3311 species are both accepted and synonyms.


# Now I need to find the accepted name for all the tips in the tree which are not accepted.
# I can do this by using the accepted_plant_name_id in the wcvp file.
# I can then use the accepted_plant_name_id to find the accepted name in the wcvp file.

# I need to find the accepted_plant_name_id for all the tips in the tree which are not accepted.
new_matching_tips_not_accepted_id <- wcvp$accepted_plant_name_id[which(new_matching_tips_not_accepted %in% wcvp$taxon_name)]
length(new_matching_tips_not_accepted_id) # 11044
length(which(is.na(new_matching_tips_not_accepted_id))) # 24 sp are missing an accepted_name_id

# So i need to drop these 24 entries from the vector and the previous vector.
cat("We have to drop these species as they have no accepted_plant_name_id:", new_matching_tips_not_accepted[which(is.na(new_matching_tips_not_accepted_id))] ,"\n")
tree <- drop.tip(tree, new_matching_tips_not_accepted[which(is.na(new_matching_tips_not_accepted_id))])
new_matching_tips_not_accepted <- new_matching_tips_not_accepted[which(!is.na(new_matching_tips_not_accepted_id))]
new_matching_tips_not_accepted_id <- new_matching_tips_not_accepted_id[which(!is.na(new_matching_tips_not_accepted_id))]


new_matching_tips_not_accepted_accepted_name <- wcvp$taxon_name[which(new_matching_tips_not_accepted_id %in% wcvp$plant_name_id)]
length(new_matching_tips_not_accepted_accepted_name) # 10843

# Are all the accepted names duplicated in wcvp?
length(new_matching_tips_not_accepted_accepted_name[duplicated(new_matching_tips_not_accepted_accepted_name)]) #21 of these are duplicated


# Finding the taxon_status of the accepted name
new_matching_tips_not_accepted_accepted_name_status <- wcvp$taxon_status[which(new_matching_tips_not_accepted_accepted_name %in% wcvp$taxon_name)]
length(which(new_matching_tips_not_accepted_accepted_name_status != "Accepted")) # 8970 of these do not point to an accepted species ??

# Creating a dataframe
rename_df_new_matches_not_accepted <- data.frame(
  taxon_name = new_matching_tips_not_accepted,
  accepted_plant_name_id = new_matching_tips_not_accepted_id,
  accepted_name = new_matching_tips_not_accepted_accepted_name,
  status_of_accept = new_matching_tips_not_accepted_accepted_name_status
)

# Okay so I need to loop through the taxon names in rename_df_new_matches_not_accepted_accepted_name and follow the accepted_plant_name_id
# to the accepted name and then check if the accepted name is accepted or not.
# I need to do this untill the named pointed to by the accepted_plant_name_id is accepted.
# for(i in seq_along(rename_df_new_matches_not_accepted$taxon_name)){ # 11020 taxon names to loop through
#   cat("Loop ", i ," :",rename_df_new_matches_not_accepted$accepted_name[i]," is ",rename_df_new_matches_not_accepted$status_of_accept[i], " \n")
#   while(rename_df_new_matches_not_accepted$status_of_accept[i] != "Accepted"){ # while the name in the accepted_name column is not accepted
#     if (is.na(rename_df_new_matches_not_accepted$accepted_plant_name_id[i]) | rename_df_new_matches_not_accepted$accepted_plant_name_id[i] == ""){ # If the accepted_plant_name_id is NA then I need to just drop the species from the tree
#       rename_df_new_matches_not_accepted$status_of_accept[i] <- "NA"
#       break
#     } else if (rename_df_new_matches_not_accepted$status_of_accept[i] == "Artificial Hybrid") {
#       rename_df_new_matches_not_accepted$status_of_accept[i] <- "Artificial Hybrid"
#       break
#     } else if (rename_df_new_matches_not_accepted$status_of_accept[i] == "Unplaced") {
#       rename_df_new_matches_not_accepted$status_of_accept[i] <- "Unplaces"
#       break
#     } else if (rename_df_new_matches_not_accepted$status_of_accept[i] == "Local Biotype") {
#       rename_df_new_matches_not_accepted$status_of_accept[i] <- "Local Biotype"
#       break
#     }
#     rename_df_new_matches_not_accepted$accepted_name[i] <- wcvp$taxon_name[which(wcvp$plant_name_id == rename_df_new_matches_not_accepted$accepted_plant_name_id[i])] # update the accepted_name column with the name pointed to by the accepted_plant_name_id
#     rename_df_new_matches_not_accepted$status_of_accept[i] <- wcvp$taxon_status[which(wcvp$plant_name_id == rename_df_new_matches_not_accepted$accepted_plant_name_id[i])] # check the status of the name pointed to by the accepted_plant_name_id
#     cat(rename_df_new_matches_not_accepted$accepted_name[i]," is ",rename_df_new_matches_not_accepted$status_of_accept[i], " \n")
#   }
# }
# saveRDS(rename_df_new_matches_not_accepted, "rename_df_new_matches_not_accepted.rds")
rename_df_new_matches_not_accepted <- readRDS("rename_df_new_matches_not_accepted.rds")


length(rename_df_new_matches_not_accepted[which(rename_df_new_matches_not_accepted$status_of_accept != "Accepted"),1]) # 149 names cannot be solved to a name which has the status accepted

# These tips need to be pruned from the dataset and from the tree.
cat("We have to drop ", length(rename_df_new_matches_not_accepted[which(rename_df_new_matches_not_accepted$status_of_accept != "Accepted"),1]),  " species as they cannot be solved to an accepted name \n")
tree <- drop.tip(tree, rename_df_new_matches_not_accepted$taxon_name[which(rename_df_new_matches_not_accepted$status_of_accept != "Accepted")])
rename_df_new_matches_not_accepted <- rename_df_new_matches_not_accepted[which(rename_df_new_matches_not_accepted$status_of_accept == "Accepted"),]


#tree$tip.label <- ifelse(tree$tip.label %in% rename_df_new_matches_not_accepted$taxon_name, rename_df_new_matches_not_accepted$accepted_name, tree$tip.label)
tree$tip.label[which(tree$tip.label %in% rename_df_new_matches_not_accepted$taxon_name)] <- rename_df_new_matches_not_accepted$accepted_name # This works because both of them is ordered


# Are all the tips in the tree now found in the WCVP_accepted file?
cat("Are all the ", length(tree$tip.label)," tips in the tree now found in the WCVP_accepted file? ", all(tree$tip.label %in% wcvp_accepted$taxon_name), "\n")

wcvp_SmB_all_accepted <- wcvp[which(wcvp$taxon_name %in% tree$tip.label),] # 10843 t
length(unique(wcvp_SmB_all_accepted$genus))
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
# So now all the tip.names in wcvp actually have the status accepted.
# Now we just need to figure out how to make all the tips in the tree have the taxon_rank == Species

# Find all the names in the tree which are found in the wcvp but which are something else than "Species"
not_species <- tree$tip.label[which(tree$tip.label %in% wcvp$taxon_name[wcvp$taxon_rank != "Species"])] # 4227 names are not species
length(not_species) 

# Fixing species that are not accepted and not species by finding their accepted name in the wcvp file
cat("There are ", length(not_species), " tips in the tree which are not species \n") # 4227


# Fixing species that are accepted but are not species
acc_not_species <- not_species[which(not_species %in% wcvp$taxon_name[wcvp$taxon_status == "Accepted"])]
length(acc_not_species) # 4227
# The accepted not species I should be able to find the species by just taking the first 2 elements of the name
tips_to_rename <- character(0)
name_for_tips_to_be_renamed <- character(0)
tips_to_drop <- character(0)
for(i in seq_along(acc_not_species)){
  tip_elements <- strsplit(acc_not_species[i], " ")[[1]]
  #print(tip_elements)
  #print(length(tip_elements))
  if(length(tip_elements) <= 2){
    tips_to_drop <- c(tips_to_drop, acc_not_species[i])
    #cat(acc_not_species[i]," Tip elements are less than 2 \n")
    next
  } else if ( any(tip_elements[1:2] == "x" | tip_elements[1:2] == "X" | tip_elements[1:2] == "×") ){
    tips_to_drop <- c(tips_to_drop, acc_not_species[i])
    #cat(acc_not_species[i]," Tip elements contains an x or an X \n")
    next
  } else if( any(tip_elements[1:2] == "Sp." | tip_elements[1:2] == "sp.")){
    tips_to_drop <- c(tips_to_drop, acc_not_species[i])
    #cat(acc_not_species[i]," Tip elements contains an Sp. or an sp. \n")
    next
  } else if ( any(is.na(tip_elements[1:2]))){
    tips_to_drop <- c(tips_to_drop, acc_not_species[i])
    #cat(acc_not_species[i]," Tip elements contains an NA \n")
    next
  }
  tips_to_rename <- c(tips_to_rename, acc_not_species[i])
  name_for_tips_to_be_renamed <- c(name_for_tips_to_be_renamed, paste(tip_elements[1:2], collapse = " "))
}

tips_to_drop
tips_to_rename

# Are all the tips to rename also accepted names?
all(name_for_tips_to_be_renamed %in% wcvp_accepted$taxon_name)
all(tips_to_drop %in% tree$tip.label)
all(tips_to_rename %in% tree$tip.label)

# Now i need to drop tips which are not tips which cannot be made into a Species name
tree <- drop.tip(tree, tips_to_drop)

# Then I need to rename all the tips which are subspecies and others into their species name
tree$tip.label[which(tree$tip.label %in% tips_to_rename)] <- name_for_tips_to_be_renamed # This works because both of them is ordered

cat("Are all the tips in the tree now both taxon_status == Accepted and taxon_rank == Species?",all(tree_test$tip.label %in% wcvp_accepted$taxon_name[which(wcvp_accepted$taxon_rank == "Species")])) # TRUE

############################################################################################################################################################################
############################################################################################################################################################################
# Now I need to take care of all the duplicate tips in the tree.
total_dup_names <- length(unique(tree$tip.label[duplicated(tree$tip.label)])) # 4964 tip names are duplicated
total_dup_names
total_dups <- length(tree$tip.label[duplicated(tree$tip.label)]) # 8587 tips in the tree have a duplicate
total_dups


# Find unique duplicated tips
unique_dup_tips <- unique(tree$tip.label[duplicated(tree$tip.label)])

# Initialize a vector to store the count of duplicates for each tip
dup_counts <- numeric(length(unique_dup_tips))

# Loop over unique duplicated tips
for (i in seq_along(unique_dup_tips)) {
  # Count occurrences in the original tree$tip.label
  dup_counts[i] <- sum(tree$tip.label == unique_dup_tips[i])
}

# Create a data frame to store results
dup_counts_df <- data.frame(tip_label = unique_dup_tips, duplicate_count = dup_counts)

head(dup_counts_df)
dup_levels <- sort(unique(dup_counts_df$duplicate_count))


# Count occurrences of each level of duplicate_count
dup_levels_counts <- table(dup_counts_df$duplicate_count)

# Create a data frame to store the results
dup_levels_counts_df <- data.frame(duplicate_count = as.numeric(names(dup_levels_counts)),
                                   count = as.numeric(dup_levels_counts))

# Print the result
print(dup_levels_counts_df)

# Create a data frame to store results
dup_levels_counts_df <- data.frame(duplicate_count = as.numeric(names(dup_levels_counts)),
                                   count = as.numeric(dup_levels_counts),
                                   Number_solved = integer(length(dup_levels_counts)),
                                   Number_not_solved = integer(length(dup_levels_counts)))



no_mono_dropped_tips <- character(0)
mono_dups <- character(0)

##############################################################################
# This function for checking duplicates is more sophisticated
# It looks trough all the duplicate tips.
# It looks for a monophyletic clade containing more than half of the number of duplicate tips.
# If it finds this clade, it drops all the tips outside the clade and drops all but one of the tips in the clade.
# If it does NOT find this clade, it drops all the duplicate tips.

# In the end it writes a small dataframe showing the number of tip names for each number of duplicates, how many of them it solved, how many of them it was unable to solve and the proportion solved.

# Current status is that I think this function runs but it takes a REALLY long time to run.
# It has to be ran on the server, which is something I havent implemented yet.

# tree_test <- tree

# if (length(tree_test$tip.label[duplicated(tree_test$tip.label)]) == 0) {
#   cat("No duplicate species names found in the tree\n")
# } else {
#   # Loop through the duplicated species names
#   list_of_dups <- unique(tree_test$tip.label[duplicated(tree_test$tip.label)])

#   for (i in seq_along(list_of_dups)) {
#     dupli_tips <- which(tree_test$tip.label == list_of_dups[i])
#     if (length(dupli_tips) > 1) {

#       # Determine half_count based on whether the number of duplicate tips is even or odd
#       half_count <- if (length(dupli_tips) %% 2 == 0) {
#         length(dupli_tips) / 2 + 1
#       } else {
#         ceiling(length(dupli_tips) / 2)
#       }

#       # Generate all combinations of duplicated tips
#       all_combinations <- lapply(half_count:length(dupli_tips), function(x) combn(dupli_tips, x, simplify = FALSE))

#       # Check if any combination forms a monophyletic clade with more than half of the duplicated tips
#       clade_formed <- any(sapply(all_combinations, function(combo) {
#         half_count <- ceiling(length(combo) / 2)
#         if (length(dupli_tips) %% 2 == 0) {
#           # For an even number of duplicate tips, require a clade larger than half the number
#           half_count <- half_count + 1
#         }
#         clade_tips <- combo[1:half_count]

#         all(sapply(clade_tips, function(tip) is.monophyletic(tree_test, tip)))
#       }))
      
#       if (clade_formed) {
#         cat("More than half of the duplicate tips of", list_of_dups[i], "form a monophyletic clade.\n")
        
#         # Keep only one tip from the first combination and remove the rest
#         if (length(all_combinations[[1]]) > 1) {
#             tree_test <- drop.tip(tree_test, tree_test$tip.label[[all_combinations[[1]][-1]]])
#         } else {
#           cat("No tips to drop from the clade.\n")
#         }

  
#         # Remove duplicate tips outside the first combination
#         non_clade_tips <- setdiff(dupli_tips, all_combinations[[1]])
#         tree_test <- drop.tip(tree_test, tree_test$tip.label[non_clade_tips])
#         mono_dups <- c(mono_dups, list_of_dups[i])
        
#         cat("Keeping only one tip from the clade and removing the rest.\n")

#         # Update the No_solved count in dup_levels_counts_df
#         idx <- which(dup_levels_counts_df$duplicate_count == length(dupli_tips))
#         dup_levels_counts_df$Number_solved[idx] <- dup_levels_counts_df$Number_solved[idx] + 1
#       } else {
#         cat("Less than half of the duplicate tips of", list_of_dups[i], "form a monophyletic clade.\n")
#         cat("Removing all duplicate tips from the tree.\n")
#         no_mono_dropped_tips <- c(no_mono_dropped_tips, list_of_dups[i])
        
#         # Update the No_not_solved count in dup_levels_counts_df
#         idx <- which(dup_levels_counts_df$duplicate_count == length(dupli_tips))
#         dup_levels_counts_df$Number_not_solved[idx] <- dup_levels_counts_df$Number_not_solved[idx] + 1
        
#         # Remove all duplicate tips from the tree
#         tree_test <- drop.tip(tree_test, tree_test$tip.label[dupli_tips])
#       }
#     }
#   }
# }



# # Print the final result
# print(dup_levels_counts_df)


# ###############################################################################
# Simple duplicate remover where I remove duplicate tips if there are more than 3
# If there is 2 I check if they are monophyletic and only remove 1 of them if they are and both if they are not
# If there is 3 duplicates I check if  any combination of the tips can form a monophyletic clade, if they can, I keep one of the species in the clade and drop the rest, if they cannot I drop all of them.

# if no duplicates are found, report it and continue
if (length(tree$tip.label[duplicated(tree$tip.label)]) == 0) {
  cat("No duplicate species names found in the tree\n")
} else {
  # Loop through the duplicated species names
  list_of_dups <- unique(tree$tip.label[duplicated(tree$tip.label)])
  for (i in seq_along(list_of_dups)) {
    cat(i,"  ",list_of_dups[i], "\n")
    dupli_tips <- which(tree$tip.label == list_of_dups[i])
    cat("Number of duplicate tips found is: ", length(dupli_tips), "\n")

    if (length(dupli_tips) == 0) {
      cat("No duplicate tips found for ", list_of_dups[i], "\n")
      next
    }

    if (length(dupli_tips) == 2) {
      if (any(dupli_tips[1] == getSisters(tree, dupli_tips[2]))) {

        cat("Duplicate species of ", list_of_dups[i], " are located next to each other in the tree\n")
        tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]])
        mono_dups <- c(mono_dups, list_of_dups[i])
        next
      } else {

        cat("Duplicate species are not located next to each other in the tree\n")
        #cat("Removing both of the species names from the tree \n")
        # I want to save the species which I had to prune due to non monophyly somewhere.
        no_mono_dropped_tips <- c(no_mono_dropped_tips, list_of_dups[i])
        tips_to_drop <- tree$tip.label[dupli_tips]
        cat("Dropping: ", tips_to_drop, "\n")
        tree <- drop.tip(tree, tips_to_drop)
        next
      }
    }

    if ( length(dupli_tips) == 3 ){
      print(dupli_tips)
      # If any 2 of the duplicated tips are located next to each other in the tree or if all 3 form a monophyletic clade then we are happy.
      # if any 2 form a monophyletic clade then we can remove the third rogue tip and remove one of the species in the mono clade.
        if (any(dupli_tips[1] == getSisters(tree, dupli_tips[2] ))) {
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[3]] )
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]] )
          mono_dups <- c(mono_dups, list_of_dups[i])
          next 

        } else if ( any(dupli_tips[1] == getSisters(tree, dupli_tips[3]))) {
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[2]] )
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]] )
          mono_dups <- c(mono_dups, list_of_dups[i])
          next

        }else if ( any(dupli_tips[2] == getSisters(tree, dupli_tips[3] ))){
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]] )
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[2]] )
          mono_dups <- c(mono_dups, list_of_dups[i])
          next
        } else {
      cat("None of the 3 tips of ", list_of_dups[i]  ,"form a monophyletic clade \n")
      no_mono_dropped_tips <- c(no_mono_dropped_tips, list_of_dups[i])
      tips_to_drop <- tree$tip.label[dupli_tips]
      cat("Dropping: ", tips_to_drop, "\n")
      tree <- drop.tip(tree, tips_to_drop)
      next
      }
    }
    if (length(dupli_tips) > 3) {
      cat("Too many duplicated tips of ",list_of_dups[i] , " removing all of them \n")
      no_mono_dropped_tips <- c(no_mono_dropped_tips, list_of_dups[i])
      tips_to_drop <- tree$tip.label[dupli_tips]
      cat("Dropping: ", tips_to_drop, "\n")
      tree <- drop.tip(tree, tips_to_drop)
      next
      }
    }
  }

 


cat("Proportion of duplicated tips solved is ", length(mono_dups)/total_dup_names, "\n")
cat("Proportion of duplicated tips not solved is ", length(no_mono_dropped_tips)/total_dup_names, "\n")
cat("Duplicated species accounted for: " , (length(mono_dups) + length(no_mono_dropped_tips))/total_dup_names, "\n")
cat("The number of duplicated tips left in the tree is: ", length(tree$tip.label[duplicated(tree$tip.label)]), "\n")



#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
wcvp_accepted_species <- wcvp_accepted[which(wcvp_accepted$taxon_rank == "Species"),]

# Are all the tip labels in the tree in the WCVP_accepted_species file?
all(tree$tip.label %in% wcvp_accepted_species$taxon_name)


# Check if all the tips in the tree are currently in the WCVP and are Accepted species names
if (length(tree$tip.label[!(tree$tip.label %in% wcvp$taxon_name[wcvp$taxon_status == "Accepted" & wcvp$taxon_rank == "Species"])]) == 0) {
  cat("All tips in the tree are in the WCVP file and are Accepted species names\n")
} else {
  cat("There are still tips in the tree that are not in the WCVP file or are not Accepted species names\n")
  cat("The tips are ", tree$tip.label[!(tree$tip.label %in% wcvp$taxon_name[wcvp$taxon_status == "Accepted" & wcvp$taxon_rank == "Species"])], "\n")
  stop("Stopping the program\n")
  break
}

# How many tips are left in the tree
cat("There are ", length(tree$tip.label), "tips left in the tree") # 64256 tips left in the tree
cat("In total we have lost ", length(no_mono_dropped_tips) + length(split_not_matchable_tips) + length(split_multi_match) + length(new_matching_tips_not_accepted), " tips from the tree")
cat(length(no_mono_dropped_tips)," due to duplicate tips not being monophyletic") 
cat(length(split_not_matchable_tips) + length(split_multi_match), " due to not being able to match the tip to the WCVP file")
cat(length(new_matching_tips_not_accepted), " due to not being able to find an accepted name for the tip")

# I should now be able to save the tree as a newick file and use it for the next step in the workflow.
cat("Saving the tree as a newick file \n")
write.tree(tree, output_file_tree)

# I also want to save the matchable tips, the split_matchable_tips and the not_matchable_tips as a data frame. Just so I can look through them when the script is done.
df_matchable_tips <- data.frame(tip = matchable_tips, match = match_name) #
df_split_matchable_tips <- data.frame(tip = split_matchable_tips, match = split_match_name) #
df_not_matchable_tips <- data.frame(tip = split_not_matchable_tips)

# Writing the data frames to files
write.table(df_matchable_tips, "matchable_tips.txt", sep = "\t", row.names = FALSE) #nolint
write.table(df_split_matchable_tips, "split_matchable_tips.txt", sep = "\t", row.names = FALSE) #nolint
write.table(df_not_matchable_tips, "not_matchable_tips.txt", sep = "\t", row.names = FALSE) #nolints


wcvp_SmB <- wcvp[which(wcvp$taxon_name %in% tree$tip.label),]
length(unique(wcvp_SmB$genus))

all(tree$tip.label %in% wcvp_accepted_species$taxon_name)

