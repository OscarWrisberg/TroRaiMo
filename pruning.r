# Load the required packages
library(ape)
library(data.table)
library(dplyr)
library(phytools)
library(stringdist)

###########################################################################################################################
# Settings for running this script locally.
#setwd("/home/au543206/GenomeDK/Trf_models/data") # Set the working directory when local
# setwd("/home/owrisberg/Trf_models/data") # Set working directory when remove
# wcvp <- readRDS("../workflow/02_adding_orders/wcvp_names_apg_aligned.rds")  # Read the WCVP names file into a data frame
# tree <- read.tree("GBMB.tre") # Read the GBMB tree
# output_file_tree <- "GBMB_pruned.tre" # Set the name of the output file
###########################################################################################################################

########################################
### Running the script on the cluster ##
########################################

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

tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

# Get tip names from the tree
tip_names <- tree$tip.label

# Find matching and non-matching tips
matching_tips <- tip_names[tip_names %in% wcvp$taxon_name] # 76935 tips are matching
not_matching_tips <- tip_names[!(tip_names %in% wcvp$taxon_name)] # Only 2939 tips are not matching

# I also need to do something about the tips that are not in the WCVP file but are in the tree.
# In an ideal world I would like to use the taxonomy_matcher to find the correct name for all the tips in the tree
# But that is for a different day.....

# Finding length of matching tips
cat("Length of matching tips ", length(matching_tips), "\n")
cat("Length of non-matching tips ", length(not_matching_tips), "\n")
# Create a data frame with tip names and families

# Hardcoding in to drop the stupid Passiflora hybrid cultivar which fuzzy matches to Passiflora hybrida.
tree <- drop.tip(tree, "Passiflora hybrid cultivar")

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
#     if (length(matches_with_distance_one) > 0) {
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

# Renaming tips in the tree with the matched tip names
tree$tip.label[which(tree$tip.label %in% matchable_tips)] <- match_name # This works because both of them is ordered


###############################################################################################################################################
# Can I loop through the unmatched tips and find the tips where we have a subsp or a variety? and where we can find the parent species.
# looping through the not_matchable_tips which when split by " " gives me a vector longer than 2.
###############################################################################################################################################

# Initialize empty vectors
split_not_matchable_tips <- character(0)
split_matchable_tips <- character(0)
split_match_name <- character(0)
split_multi_match <- character(0)

# Loop through the not_matchable_tips
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

# Are all species accounted for?
cat("Are all species accounted for: "all(not_matchable_tips %in% split_not_matchable_tips | not_matchable_tips %in% split_matchable_tips| not_matchable_tips %in% split_multi_match))

Saving vectors as RDS files for easy loading if I need to rerun the script
saveRDS(split_not_matchable_tips, "split_not_matchable_tips.rds")
saveRDS(split_matchable_tips, "split_matchable_tips.rds")
saveRDS(split_match_name, "split_match_name.rds")
saveRDS(split_multi_match, "split_multi_match.rds")

# Loading RDS files
# split_not_matchable_tips <- readRDS("split_not_matchable_tips.rds")
# split_matchable_tips <- readRDS("split_matchable_tips.rds")
# split_match_name <- readRDS("split_match_name.rds")
# split_multi_match <- readRDS("split_multi_match.rds")

# And now we can again rename the tips based on the matches we found
cat("Renaming tips we found during split names approach \n")
tree$tip.label[which(tree$tip.label %in% split_matchable_tips)] <- split_match_name

# Remove tips that are still not matched
cat("Removing tips that are still not matched \n")
tree <- drop.tip(tree, split_not_matchable_tips)

tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

no_mono_dropped_tips <- character(0)
total_dups <- length(tree$tip.label[duplicated(tree$tip.label)])
mono_dups <- character(0)

# if no duplicates are found, report it and continue
if (length(tree$tip.label[duplicated(tree$tip.label)]) == 0) {
  cat("No duplicate species names found in the tree\n")
} else {
  # Loop through the duplicated species names
  for (i in seq_along(tree$tip.label[duplicated(tree$tip.label)])) {
    #print(tree$tip.label[duplicated(tree$tip.label)][i])
    dupli_tips <- which(tree$tip.label == tree$tip.label[duplicated(tree$tip.label)][i])
    #cat("Number of duplicate tips found is: ", length(dupli_tips), "\n")

    if (length(dupli_tips) == 2) {
      if (dupli_tips[1] == getSisters(tree, dupli_tips[2])) {

        cat("Duplicate species of ", tree$tip.label[duplicated(tree$tip.label)][i], " are located next to each other in the tree\n")
        cat("Removing one of the species names from the tree \n")
        tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]])
        mono_dups <- c(mono_dups, tree$tip.label[duplicated(tree$tip.label)][i])
        next
      } else {

        #cat("Duplicate species are not located next to each other in the tree\n")
        #cat("Removing both of the species names from the tree \n")
        # I want to save the species which I had to prune due to non monophyly somewhere.
        no_mono_dropped_tips <- c(no_mono_dropped_tips, tree$tip.label[duplicated(tree$tip.label)][i])
        tree <- drop.tip(tree, tree$tip.label[duplicated(tree$tip.label)][i])
        next
      }
    }

    if ( length(dupli_tips == 3 )){
      # If any 2 of the duplicated tips are located next to each other in the tree or if all 3 form a monophyletic clade then we are happy.
      # if any 2 form a monophyletic clade then we can remove the third rogue tip and remove one of the species in the mono clade.
        if (dupli_tips[1] == getSisters(tree, dupli_tips[2] )) {
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[3]] )
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]] )
          mono_dups <- c(mono_dups, tree$tip.label[duplicated(tree$tip.label)][i])
          next 

        } else if ( dupli_tips[1] == getSisters(tree, dupli_tips[3])) {
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[2]] )
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]] )
          mono_dups <- c(mono_dups, tree$tip.label[duplicated(tree$tip.label)][i])
          next

        }else if ( dupli_tips[2] == getSisters(tree, dupli_tips[3] )){
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[1]] )
          tree <- drop.tip(tree, tree$tip.label[dupli_tips[2]] )
          mono_dups <- c(mono_dups, tree$tip.label[duplicated(tree$tip.label)][i])
          next
        }
    } else {
      cat("None of the 3 tips of ", tree$tip.label[duplicated(tree$tip.label)][i]  ,"form a monophyletic clade \n")
      no_mono_dropped_tips <- c(no_mono_dropped_tips, tree$tip.label[duplicated(tree$tip.label)][i])
      tree <- drop.tip(tree, tree$tip.label[duplicated(tree$tip.label)][i])
      next
    }

    if (length(dupli_tips > 3)) {
      cat("Too many duplicated tips of ",tree$tip.label[duplicated(tree$tip.label)][i] , " removing all of them \n")
      no_mono_dropped_tips <- c(no_mono_dropped_tips, tree$tip.label[duplicated(tree$tip.label)][i])
      tree <- drop.tip(tree, tree$tip.label[duplicated(tree$tip.label)][i])
      next
    }
  }
}

cat("Proportion of duplicated tips solved is ", length(mono_dups)/total_dups, "\n")
cat("Proportion of duplicated tips not solved is ", length(no_mono_dropped_tips)/total_dups, "\n")
cat("Duplicated species accounted for: " , length(mono_dups) + length(no_mono_dropped_tips), "\n")


# Check if there are any tips left in the tree that are not in the WCVP file
if (length(tree$tip.label[!(tree$tip.label %in% wcvp$taxon_name)]) == 0) {
  cat("All tips in the tree are in the WCVP file\n")
} else {
  cat("There are still tips in the tree that are not in the WCVP file\n")
  cat("The tips are ", tree$tip.label[!(tree$tip.label %in% wcvp$taxon_name)], "\n")
  stop("Stopping the program\n")
  break
}

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
