##########################################################################################################################
########################################### Load the required packages ###################################################
##########################################################################################################################

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "castor", "dplyr", "stringdist")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


###########################################################################################################################
######################################## Settings for running this script locally. ########################################
###########################################################################################################################

# setwd("/home/au543206/GenomeDK/Trf_models/data") # Set the working directory when local
# wcvp <- readRDS("../workflow/02_adding_orders/wcvp_names_apg_aligned.rds")  # Read the WCVP names file into a data frame
# tree <- read.tree("GBMB.tre") # Read the GBMB tree
# output_file_tree <- "GBMB_pruned.tre" # Set the name of the output file

###########################################################################################################################
########################################### Running the script on the cluster #############################################
###########################################################################################################################

# # Command line arguments
input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
input_file_wcvp <- commandArgs(trailingOnly = TRUE)[2]
output_file_tree <- commandArgs(trailingOnly = TRUE)[3]

# Read the WCVP names file into a data frame
cat("Opening ", input_file_wcvp, "\n")  
wcvp <- readRDS(input_file_wcvp)

# Read the GBMB tree
cat("Opening ", input_file_tree, "\n")
tree <- read.tree(input_file_tree)

###########################################################################################################################
########################################################### Basic Setup  ##################################################
###########################################################################################################################

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
cat("Length of matching tips ", length(matching_tips), "\n") #76935
cat("Length of non-matching tips ", length(not_matching_tips), "\n") #2938
cat("Length of matching tips accepted ", length(matching_tips_accepted), "\n") #66166
cat("Length of matching tips not accepted ", length(matching_tips_not_accepted), "\n") #14355

# Create a data frame with tip names and families

# Looking at the not_matching_tips I need to loop through them in order to check if they are indeed unmatchable or if they can be matched.
# I could use the taxonomy_matcher to find the correct name for all the tips in the tree, but I think this would be too time consuming.
# Instead I am just going to loop through the names and see if I can find a match in the WCVP file using the grep function.

# Initialize empty vectors
not_matchable_tips <- character(0)
matchable_tips <- character(0)
match_name <- character(0)

# Loop through the not_matching_tips
for (i in seq_along(not_matching_tips)) {
  # Writing a progress bar
  if (!i %% 50) cat("Percentage done", format(round((i / length(not_matching_tips)) * 100, 2), nsmall = 2), " at ", format(Sys.time(), '%H:%M:%S'), "\n")

  # Check for exact match
  if (not_matching_tips[i] %in% wcvp$taxon_name) {
    #cat("Exact match found for ", not_matching_tips[i], "\n")
    matchable_tips <- c(matchable_tips, not_matching_tips[i])
    match_name <- c(match_name, not_matching_tips[i])
  } else {
    # Check for matches with one substitution, insertion, or deletion
    potential_matches <- stringdist::stringdistmatrix(not_matching_tips[i], wcvp$taxon_name, method = "lv", useNames = TRUE)
    
    # Filter for matches with distance 1
    matches_with_distance_one <- wcvp$taxon_name[potential_matches <= 1]
    
    # Check if any matches were found
    if (length(matches_with_distance_one) > 0 & length(matches_with_distance_one) == 1) {
      #cat("Match found for ", not_matching_tips[i], " with one substitution, insertion, or deletion\n")
	    if (length(matches_with_distance_one) == 1) {
        #cat("Match is ", matches_with_distance_one[1], "\n")
        matchable_tips <- c(matchable_tips, not_matching_tips[i])
        match_name <- c(match_name, matches_with_distance_one[1])
	    } else {
		    cat("Multiple matches found for ", not_matching_tips[i], "\n")
		    cat("Matches are ", matches_with_distance_one, "\n")
		    not_matchable_tips <- c(not_matchable_tips, not_matching_tips[i])
      } 
	  } else {
      # No match found
      #cat("No match found for ", not_matching_tips[i], "\n")
      not_matchable_tips <- c(not_matchable_tips, not_matching_tips[i])
    }
  }
}

# Saving vectors as RDS files
saveRDS(not_matchable_tips, "not_matchable_tips_1.rds") 
saveRDS(matchable_tips, "matchable_tips_1.rds") 
saveRDS(match_name, "match_name_1.rds")

# # Loading RDS files
not_matchable_tips <- readRDS("not_matchable_tips_1.rds")
matchable_tips <- readRDS("matchable_tips_1.rds")
match_name <- readRDS("match_name_1.rds")


cat("Are all the matchable tips found in the tree? ", all(matchable_tips %in% tree$tip.label), " and are all the match_names found in the wcvp$taxon_name", all(match_name %in% wcvp$taxon_name), "\n")


# Find the indices of tips in the tree that match the split_matchable_tips
matching_indices <- match(tree$tip.label, matchable_tips)

# Replace the matched tips with split_match_name
tree$tip.label[!is.na(matching_indices)] <- match_name

# Are all species accounted for?
cat("Have all the match_name's been incorporated into the tip.labels: ",all(match_name %in% tree$tip.label))
cat("Are all species accounted for: ", all(not_matching_tips %in% not_matchable_tips | not_matching_tips %in% matchable_tips), "\n")

###############################################################################################################################################
# Can I loop through the unmatched tips and find the tips where we have a subsp or a variety? and where we can find the parent species.
# looping through the not_matchable_tips which when split by " " gives me a vector longer than 2.
###############################################################################################################################################

# Initialize empty vectors
split_not_matchable_tips <- character(0)
split_matchable_tips <- character(0)
split_match_name <- character(0)
split_multi_match <- character(0)

#Loop through the not_matchable_tips
for (i in seq_along(not_matchable_tips)) {
  # Writing a progress bar
  if (!i %% 50) cat("Percentage done", format(round((i / length(not_matchable_tips)) * 100, 2), nsmall = 2), " at ", format(Sys.time(), '%H:%M:%S'), "\n")

  # Split the tip by space
  tip_elements <- strsplit(not_matchable_tips[i], " ")[[1]]
  
  # If the split list is only of length 2, then we can go to next tip
  if(length(tip_elements) <= 2){
    split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
	next
  }

  if( any(tip_elements == "x" | tip_elements == "X") ){
    # If the tip contains an x or an X then we can go to next tip
    #cat("Cannot match ", not_matchable_tips[i], " because it contains an x or an X \n")
    split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
    next
  } 

  if ( any(tip_elements == "Sp." | tip_elements == "sp.")){
    # If the tip contains an Sp. or an sp. then we can go to next tip
    #cat("Cannot match ", not_matchable_tips[i], " because it contains an Sp. or an sp. \n")
    split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
    next
  }

  # Otherwise we can try to find a match for the first two elements
  if (length(tip_elements) > 2) {
    # Take the first two elements
    tip_to_search <- paste(tip_elements[1:2], collapse = " ")
    
    # Check for exact match
    if (tip_to_search %in% wcvp$taxon_name) {
      #cat("Exact match found for ", not_matchable_tips[i], "\n")
      split_matchable_tips <- c(split_matchable_tips, not_matchable_tips[i])
      split_match_name <- c(split_match_name, tip_to_search)
    } else {
      # Check for matches with one substitution, insertion, or deletion
      potential_matches <- stringdist::stringdistmatrix(tip_to_search, wcvp$taxon_name, method = "lv", useNames = TRUE)
      
      # Filter for matches with distance 1
      matches_with_distance_one <- wcvp$taxon_name[potential_matches <= 1]
      
      # Check if any matches were found
      if (length(matches_with_distance_one) > 0 & length(matches_with_distance_one) == 1) {
        #cat("Match found for ", not_matchable_tips[i], " with one substitution, insertion, or deletion\n")
        #cat("Match is ", matches_with_distance_one[1], "\n")
        split_matchable_tips <- c(split_matchable_tips, not_matchable_tips[i])
        split_match_name <- c(split_match_name, matches_with_distance_one[1])
      }

       if(length(matches_with_distance_one) > 0 & length(matches_with_distance_one) > 1) {
        #cat("Multiple matches found for ", not_matchable_tips[i], "\n")
        #cat("Matches are ", matches_with_distance_one, "\n")
        split_multi_match <- c(split_multi_match, not_matchable_tips[i])

      } else {
        # No match found
        #cat("No match found for ", not_matchable_tips[i], "\n")
        split_not_matchable_tips <- c(split_not_matchable_tips, not_matchable_tips[i])
      }
    }
  }
}

# Saving vectors as RDS files for easy loading if I need to rerun the script
cat("Saving Rds files \n")
saveRDS(split_not_matchable_tips, "split_not_matchable_tips.rds")
saveRDS(split_matchable_tips, "split_matchable_tips.rds")
saveRDS(split_match_name, "split_match_name.rds")
saveRDS(split_multi_match, "split_multi_match.rds")
cat("Rds files saved \n")

## Loading RDS files
split_not_matchable_tips <- readRDS("split_not_matchable_tips.rds")
split_matchable_tips <- readRDS("split_matchable_tips.rds")
split_match_name <- readRDS("split_match_name.rds")
split_multi_match <- readRDS("split_multi_match.rds")


# Are all species accounted for?
cat("Are all species accounted for: ", all(not_matchable_tips %in% split_not_matchable_tips | not_matchable_tips %in% split_matchable_tips| not_matchable_tips %in% split_multi_match), "\n")

# And now we can again rename the tips based on the matches we found
cat("Number of tips in split_matchable_tips is: ", length(split_matchable_tips), "\n")
cat("Number of tips in split_match_name is: ", length(split_match_name), "\n")
cat("Number of tips in tree which are found in matchable names are: ", length(tree$tip.label[which(tree$tip.label %in% split_matchable_tips)]))

# Removing the row from split_matchable_tips and split_match_name if the tip is not in the tree
split_match_name <- split_match_name[which(split_matchable_tips %in% tree$tip.label)]
split_matchable_tips <- split_matchable_tips[which(split_matchable_tips %in% tree$tip.label)]


# Find the indices of tips in the tree that match the split_matchable_tips
matching_indices <- match(tree$tip.label, split_matchable_tips)

# Replace the matched tips with split_match_name
tree$tip.label[!is.na(matching_indices)] <- split_match_name

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
matching_tips <- tip_names[tip_names %in% wcvp$taxon_name] # 78208 tips are matching
not_matching_tips <- tip_names[!(tip_names %in% wcvp$taxon_name)] # Only 0 tips are not matching
matching_tips_accepted <- tip_names[tip_names %in% wcvp$taxon_name[wcvp$taxon_status == "Accepted"]] # 67164 out of the matching tips are accepted
matching_tips_not_accepted <- tip_names[tip_names %in% wcvp$taxon_name[wcvp$taxon_status != "Accepted"]] # 14355 out of the matching tips are not accepted


# Finding length of matching tips
cat("Length of matching tips ", length(matching_tips), "\n") # 78208 we have gained 1273 new matching species
cat("Length of non-matching tips ", length(not_matching_tips), "\n") #0
cat("Length of matching tips accepted ", length(matching_tips_accepted), "\n") #67164 we have gained 998 new accepted species
cat("Length of matching tips not accepted ", length(matching_tips_not_accepted), "\n") #14673 we have gained 318 new not accepted species

# End of dealing with tips which were not matching to the WCVP file.
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
wcvp_n_accepted <- wcvp[which(wcvp$taxon_status != "Accepted")]

# Because some taxon_names are both accepted and not accepted based on different authors, we need to find all the taxon_names which dont have a version which is accepted.
new_matching_tips_not_accepted_old <- unique(tree$tip.label[which(!tree$tip.label %in% wcvp_accepted$taxon_name)])
length(new_matching_tips_not_accepted_old) # 10984

# Stating the number of species which are both accepted and not accepted at the same time.
cat("There are ",length(matching_tips_not_accepted)-length(new_matching_tips_not_accepted_old) , " tips in the tree which are both Accepted and Synonym due to different authors. \n") # 3689
cat("This means there are ",length(new_matching_tips_not_accepted_old)," tips in the tree which are not accepted at all \n") # still 10984 tips are not accepted, which is down from 14355. This means 3689 species are both accepted and synonyms.

# Now I need to find the accepted name for all the tips in the tree which are not accepted.
# I can do this by using the accepted_plant_name_id in the wcvp file.
# I can then use the accepted_plant_name_id to find the accepted name in the wcvp file.

# Could it be because there are more than one row where wcvp$taxon_name matches new_matching_tips_not_accepted YES THIS SEEMS TO BE THE CASE!
# I could test this with a for loop.
# This for loop loops through all the names in new_matching_tips_not_accepted and checks
# if there are more than one row in the wcvp file where the taxon_name matches the name in new_matching_tips_not_accepted
# If there is more than one row, it checks if there is only one accepted_plant_name_id which is not NA
          # If there is only one accepted_plant_name_id which is not NA, then it prints the name and the accepted_plant_name_id
          # But this is not the case for any of the tips
#Otherwise it adds the tip to a vector called test_tip_labels_more_than_1_match
# Which is then dropped from the tree.

new_matching_tips_not_accepted <- unique(tree$tip.label[which(!tree$tip.label %in% wcvp_accepted$taxon_name)])
test_tip_labels_more_than_1_match <- character(0)
for(i in seq_along(new_matching_tips_not_accepted)){
  if (length(which(wcvp$taxon_name %in% new_matching_tips_not_accepted[i])) > 1) {
    #print(new_matching_tips_not_accepted[i])
    rows_with_name <- wcvp[which(wcvp$taxon_name %in% new_matching_tips_not_accepted[i])]
    rows_with_name_id <- rows_with_name$accepted_plant_name_id
    #print(rows_with_name_id)
    if(length(rows_with_name_id)-length(is.na(rows_with_name_id)) == 1) {
      cat("This one can be fixed") # None of the tips can be fixed
      print(rows_with_name_id)
      print(!is.na(rows_with_name_id))
    } else {
    #print(new_matching_tips_not_accepted[i])
    test_tip_labels_more_than_1_match <- c(test_tip_labels_more_than_1_match, new_matching_tips_not_accepted[i])
    }
  }
}
saveRDS(test_tip_labels_more_than_1_match, "test_tip_labels_more_than_1_match.rds")
test_tip_labels_more_than_1_match <- readRDS("test_tip_labels_more_than_1_match.rds")

# I guess I have to remove these tips from the tree as I am unable to find their accepted name
tree <- drop.tip(tree, test_tip_labels_more_than_1_match)

# # Updating the tips in the tree which are not accepted
# new_matching_tips_not_accepted <- unique(tree$tip.label[which(!tree$tip.label %in% wcvp_accepted$taxon_name)]) # is this unique here necessary?

# # Finding their ID's
# new_matching_tips_not_accepted_id <- wcvp$accepted_plant_name_id[which(wcvp$taxon_name %in% new_matching_tips_not_accepted)]
# length(new_matching_tips_not_accepted_id) # 11612 

# # Now I need to deal with the species where the matching id is NA
# length(which(is.na(new_matching_tips_not_accepted_id))) # 189 sp are missing an accepted_name_id

# # So i need to drop these 140 entries from the vector and the previous vector.
# cat("We have to drop these species as they have no accepted_plant_name_id: \n")
# for( i in seq_along(new_matching_tips_not_accepted[which(is.na(new_matching_tips_not_accepted_id))])) {
#   print(new_matching_tips_not_accepted[which(is.na(new_matching_tips_not_accepted_id))][i])
# }

# # First I drop the species from the tree
# tree <- drop.tip(tree, new_matching_tips_not_accepted[which(is.na(new_matching_tips_not_accepted_id))])

#new_matching_tips_not_accepted <- new_matching_tips_not_accepted_old[which(!is.na(new_matching_tips_not_accepted_id))] # This is a list of all the names.

# This finds the names and the accepted_plant_name_id for the names which are not accepted
new_matching_tips_not_accepted_name_id <- wcvp[which(wcvp$taxon_name %in% new_matching_tips_not_accepted_old),c("taxon_name", "accepted_plant_name_id", "taxon_status")] # This line of code works to find the rows where the names are present in wcvp
dim(new_matching_tips_not_accepted_name_id) # 11612-3


# This line of code works to find the rows where the names are present in wcvp
new_matching_tips_not_accepted_accepted_name <- wcvp[which(wcvp$plant_name_id %in% new_matching_tips_not_accepted_name_id$accepted_plant_name_id),c("taxon_name", "plant_name_id", "taxon_status")]
colnames(new_matching_tips_not_accepted_accepted_name)[which(colnames(new_matching_tips_not_accepted_accepted_name) == "taxon_name")] <- "accepted_name"
colnames(new_matching_tips_not_accepted_accepted_name)[which(colnames(new_matching_tips_not_accepted_accepted_name) == "taxon_status")] <- "accepted_name_taxon_status"
length(new_matching_tips_not_accepted_name_id$accepted_plant_name_id) # 11415


# Database of all the names which are not accepted and their accepted_plant_name_id
new_matching_tips_df_test <- merge(new_matching_tips_not_accepted_name_id, new_matching_tips_not_accepted_accepted_name, by.x = "accepted_plant_name_id", by.y = "plant_name_id", all.x = TRUE)
dim(new_matching_tips_df_test) # 11612 rows and 5 columns

new_matching_tips_df_test[which(new_matching_tips_df_test$accepted_name_taxon_status != "Accepted")] # 46 rows and 3 columns which we need to remove because the status is not accepted
new_matching_tips_df_test[which(is.na(new_matching_tips_df_test$accepted_plant_name_id))] # 189 which we need to remove because the accepted_plant_name_id is NA
length(tree$tip.label) # 77685

# I need to drop these 235 entries from the tree.
# Dropping tips where the accepted_plant_name_id is NA
tree <- drop.tip(tree, new_matching_tips_df_test$taxon_name[which(is.na(new_matching_tips_df_test$accepted_plant_name_id))])

# Dropping the tips where the accepted_name_taxon_status is not Accepted
tree <- drop.tip(tree, new_matching_tips_df_test$taxon_name[which(new_matching_tips_df_test$accepted_name_taxon_status != "Accepted")])

#
new_matching_tips_df_test_tree <- new_matching_tips_df_test[which(new_matching_tips_df_test$taxon_name %in% tree$tip.label),]


# Save the original tip labels
original_tip_labels <- tree$tip.label

# Create an empty data frame to store changes
changes_df <- data.frame(original_name = character(0), new_name = character(0), stringsAsFactors = FALSE)

# Loop through the tips and update the tree
for (i in seq_len(nrow(new_matching_tips_df_test_tree))) {
  original_name <- new_matching_tips_df_test_tree$taxon_name[i]
  new_name <- new_matching_tips_df_test_tree$accepted_name[i]
  
  # Find the index of the tip in the tree
  tip_index <- which(tree$tip.label == original_name)
  
  if (length(tip_index) > 0) {
    # Update the tip label in the tree
    tree$tip.label[tip_index] <- new_name
    
    # Save the change to the data frame
    changes_df <- rbind(changes_df, data.frame(original_name = original_name, new_name = new_name, stringsAsFactors = FALSE))
  }
}

# Print the changes
#print(changes_df)

# Check if tip labels were changed
changed_tip_labels <- tree$tip.label != original_tip_labels

# Show the changed tip labels
length(which(changed_tip_labels == TRUE)) # 10355 tips have had their names changed

# Are all the tips in the tree now found in the WCVP_accepted file?
cat("Are all the ", length(tree$tip.label)," tips in the tree now found in the WCVP_accepted file? ", all(tree$tip.label %in% wcvp_accepted$taxon_name), "\n") # YES I am so happy!

write.tree(tree, "GBMB_all_accepted_names.tre") # Everything here looks fine. The species are in their correct spots

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
# So now all the tip.names in wcvp actually have the status accepted.
# Now we just need to figure out how to make all the tips in the tree have the taxon_rank == Species

# Find all the names in the tree which are found in the wcvp but which are something else than "Species"
not_species <- tree$tip.label[which(tree$tip.label %in% wcvp_accepted$taxon_name[wcvp_accepted$taxon_rank != "Species"])] # 3859 names are not species
length(not_species) 

# Fixing species that are not accepted and not species by finding their accepted name in the wcvp file
cat("There are ", length(not_species), " tips in the tree which are not species \n") # 3859

# Fixing species that are accepted but are not species
acc_not_species <- not_species[which(not_species %in% wcvp_accepted$taxon_name[wcvp_accepted$taxon_status == "Accepted"])]
length(acc_not_species) # 3859

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


# Are all the tips to rename also accepted names?
all(name_for_tips_to_be_renamed %in% wcvp_accepted$taxon_name)
all(tips_to_drop %in% tree$tip.label)
all(tips_to_rename %in% tree$tip.label)
length(name_for_tips_to_be_renamed) == length(tips_to_rename)

# Now i need to drop tips which are not tips which cannot be made into a Species name
tree <- drop.tip(tree, tips_to_drop)

# Find the indices of tips in the tree that match the split_matchable_tips
matching_indices <- match(tree$tip.label, tips_to_rename)

# Replace the matched tips with split_match_name
tree$tip.label[!is.na(matching_indices)] <- name_for_tips_to_be_renamed

cat("Are all ",length(tree$tip.label)," tips in the tree now both taxon_status == Accepted and taxon_rank == Species?",all(tree$tip.label %in% wcvp_accepted$taxon_name[which(wcvp_accepted$taxon_rank == "Species")])) # TRUE

write.tree(tree, "GBMB_all_accepted_all_species.tre") # Everything here looks fine. The species are in their correct spots

############################################################################################################################################################################
############################################################################################################################################################################

# Now I need to take care of all the duplicate tips in the tree.
total_dup_names <- length(unique(tree$tip.label[duplicated(tree$tip.label)])) # 4855 tip names are duplicated
total_dups <- length(tree$tip.label[duplicated(tree$tip.label)]) # 7970 tips in the tree have a duplicate

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

print(dup_counts_df)
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

# ###############################################################################
# Simple duplicate remover where I remove duplicate tips if there are more than 3
# If there is 2 I check if they are monophyletic and only remove 1 of them if they are and both if they are not
# If there is 3 duplicates I check if  any combination of the tips can form a monophyletic clade, if they can, I keep one of the species in the clade and drop the rest, if they cannot I drop all of them.

length_before_dup_removal <- length(tree$tip.label)

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

 length_after_dup_removal <- length(tree$tip.label)

cat("Tips lost due to them being duplicates is: ", length_before_dup_removal - length_after_dup_removal, "\n") # 10137!!! My fucking god that is many!
cat("Proportion of duplicated tips solved is ", length(mono_dups)/total_dup_names, "\n")
cat("Proportion of duplicated tips not solved is ", length(no_mono_dropped_tips)/total_dup_names, "\n")
cat("The number of duplicated tips left in the tree is: ", length(tree$tip.label[duplicated(tree$tip.label)]), "\n")
cat("The number of tips in the tree is ", length(tree$tip.label))

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
cat("There are ", length(tree$tip.label), "tips left in the tree \n") # 64256 tips left in the tree
cat("In total we have lost ", length(tip_names) - length(tree$tip.label) , " tips from the tree")

# I should now be able to save the tree as a newick file and use it for the next step in the workflow.
cat("Saving the tree as a newick file \n")
write.tree(tree, output_file_tree)



#############################################################################
# Area with old and outdated code which I might want to use later


# # Okay so I need to loop through the taxon names in rename_df_new_matches_not_accepted_accepted_name and follow the accepted_plant_name_id
# # to the accepted name and then check if the accepted name is accepted or not.
# # I need to do this untill the named pointed to by the accepted_plant_name_id is accepted.

# # I will create a small dataframe to keep track of how many of the Non-Accepted names cannot be changed to an Accepted name and the reason why.
# df_finding_accepted_ids <- data.frame(
#   found_accepted_name = numeric(0),
#   found_na = numeric(0),
#   found_artificial_hybrid = numeric(0),
#   found_unplaced = numeric(0),
#   found_local_biotype = numeric(0)
# )


# # Running the while loop for each species which is not taxon_status = Accepted
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
#       rename_df_new_matches_not_accepted$status_of_accept[i] <- "Unplaced"
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





# # I also want to save the matchable tips, the split_matchable_tips and the not_matchable_tips as a data frame. Just so I can look through them when the script is done.
# df_matchable_tips <- data.frame(tip = matchable_tips, match = match_name) #
# df_split_matchable_tips <- data.frame(tip = split_matchable_tips, match = split_match_name) #
# df_not_matchable_tips <- data.frame(tip = split_not_matchable_tips)

# # Writing the data frames to files
# write.table(df_matchable_tips, "matchable_tips.txt", sep = "\t", row.names = FALSE) #nolint
# write.table(df_split_matchable_tips, "split_matchable_tips.txt", sep = "\t", row.names = FALSE) #nolint
# write.table(df_not_matchable_tips, "not_matchable_tips.txt", sep = "\t", row.names = FALSE) #nolints



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
