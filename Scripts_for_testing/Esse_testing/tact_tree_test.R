# This script is supposed to load one of the Tacted trees that Melanie made
# Then I want to cut it into subtrees based on the orders and then save each of these subtrees as a new tree
# The idea is that ClaDs is struggeling alot with trees where the sampling fraction is very low

# Because all the species should be in the Tacted trees I want to see if it is easier for handle the entire tree for each order.

# Setting order in question
order_in_question <- "Laurales"

# Load libraries
library(ape)
library(phytools)
library(data.table)

# Set working directory
setwd("/home/au543206/Documents/Tact_trees")




############################################################################################################
####################  Loading the tree and cutting it into subtrees  #######################################
############################################################################################################
# Load tree
tree <- read.tree("gbmb_matched_monophyletic_orders_63306956-323.tacted.newick.tre")

# remove '' from tip labels
tree$tip.label <- gsub("'", "", tree$tip.label)

# print 1000 random tip labels
print(tree$tip.label[sample(1:length(tree$tip.label), 1000)])

# Create a dataframe with all the tip labels
df_tips <- data.frame(tip_label = tree$tip.label)
df_tips$no_text <- gsub("[^0-9]", "", df_tips$tip_label)

# remove anything that is not numbers from the tip lables
tree$tip.label <- gsub("[^0-9]", "", tree$tip.label)

# Load wcvp
# The problem with Melanies Tacted trees is that wcvp have changed the wcvp_index since their special edition
# This means that if you download the newest issue of 
wcvp <- readRDS("/home/au543206/Documents/world_checklist/World_checklist_downloads/01_24_2024/wcvp_names_apg_aligned.rds")
wcvp_old <- fread("/home/au543206/Documents/world_checklist/checklist_names.txt", sep = "|")
wcvp_special_edition <- fread("/home/au543206/Documents/world_checklist/World_checklist_downloads/Special_edition/wcvp_names.txt", sep = "|")

# Check if all the tips are in wcvp old
df_tips$in_wcvp_old <- df_tips$tip_label %in% wcvp_old$plant_name_id

# Check if all the tips are in wcvp
df_tips$in_wcvp <- df_tips$no_text %in% wcvp$plant_name_id

# Check if all the tips are in wcvp special edition
df_tips$in_wcvp_special_edition <- df_tips$tip_label %in% wcvp_special_edition$plant_name_id

# Can we chech how many tips are represented in atleast one of the versions of wcvp
df_tips$in_wcvp_any <- df_tips$in_wcvp | df_tips$in_wcvp_old | df_tips$in_wcvp_special_edition

# Summarize how many tips are in old and how many are not
table(df_tips$in_wcvp_old)
table(df_tips$in_wcvp)
table(df_tips$in_wcvp_special_edition)
table(df_tips$in_wcvp_any)

# print the species that are not present in any of the versions of wcvp
df_tips$tip_label[which(df_tips$in_wcvp_any == FALSE)]

# If a species is found in a specific version of the wcvp then we add that versions species name to the df_tips data frame
df_tips$wcvp_old <- NA
df_tips$wcvp_old[which(df_tips$in_wcvp_old == TRUE)] <- wcvp_old$taxon_name[match(df_tips$tip_label[which(df_tips$in_wcvp_old == TRUE)], wcvp_old$plant_name_id)]

df_tips$wcvp <- NA
df_tips$wcvp[which(df_tips$in_wcvp == TRUE)] <- wcvp$taxon_name[match(df_tips$no_text[which(df_tips$in_wcvp == TRUE)], wcvp$plant_name_id)]

df_tips$wcvp_special_edition <- NA
df_tips$wcvp_special_edition[which(df_tips$in_wcvp_special_edition == TRUE)] <- wcvp_special_edition$taxon_name[match(df_tips$tip_label[which(df_tips$in_wcvp_special_edition == TRUE)], wcvp_special_edition$plant_name_id)]

# find the 19 tips which have no match in the special edition
df_tips[which(is.na(df_tips$wcvp_special_edition)),]

# The two other dataframes agree on a name for these 19 species so Ill just add them to the wcvp_special_edition dataframe column
df_tips$wcvp_special_edition[which(is.na(df_tips$wcvp_special_edition))] <- df_tips$wcvp[which(is.na(df_tips$wcvp_special_edition))]

# All the tip labels are in the wcvp_special_edition dataframe
all(tree$tip.label %in% df_tips$wcvp_special_edition)

# Rename all the tips in the tree based on the wcvp_special edition names and if they are missing use the wcvp names
tree$tip.label <- df_tips$wcvp_special_edition

all(tree$tip.label %in% wcvp$taxon_name)
length(which(tree$tip.label %in% wcvp$taxon_name))

# find the length of the tips which are not in wcvp
length(which(!tree$tip.label %in% wcvp$taxon_name)) # 536 tips are not in wcvp
tips_not_in_wcvp <- tree$tip.label[which(!tree$tip.label %in% wcvp$taxon_name)]
length(tree$tip.label) # 331505 tip in total


# print the tips which are not in wcvp
tips_not_in_wcvp

# Are these found in the special issue?
wcvp_special_edition[which(wcvp_special_edition$taxon_name %in% tips_not_in_wcvp),] # Yes they are all found in the special issue

# Now how do we 


# Can you find all the rows in wcvp which have the name Machilus in the genus column
sort(wcvp$taxon_name[which(wcvp$genus == "Machilus")])

# Filter the wcvp dataset to include only rows where taxon_status == "Accepted"
wcvp_accepted <- subset(wcvp, taxon_status == "Accepted")
wcvp_accepted_species <- subset(wcvp_accepted, taxon_rank == "Species")

# Load apg file
apg <- read.csv("/home/au543206/Documents/world_checklist/World_checklist_downloads/01_24_2024/apgweb_parsed.csv")

# Find unique families
unique_families <- unique(wcvp_accepted_species$family)
unique_families <- as.character(unique_families)

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))
non_mono_family <- character(0)

####################################################################################
####################  Finding the order for each family  ###########################
####################################################################################

# Find order function
find_order <- function(fams, apg) {
  fam_list <- character(0)
  orders <- character(0)

  for (family in unique(fams)) {
    order <- apg[which(apg$Syn_Fam == family), "Clade"] # Finding the order of that family in APG file
    order <- as.character(order[1]) # Selecting the order of the family
    cat("family ", family, "Order", order, "\n") # Printing the family and order
    fam_list <- c(fam_list, family)
    orders <- c(orders, order)
  }

  df_orders <- data.frame(family = fam_list, order = orders)
  return(df_orders)
}

# Running the function
cat("Finding the order for each of the families \n")
family_orders <- find_order(unique_families, apg)
length(family_orders$order)

# Merging the wcvp and family_orders data frames
cat("Merging the wcvp and family_orders data frames \n")
wcvp_accepted_species_orders <- merge(wcvp_accepted_species, family_orders, by.x = "family", by.y = "family")

# Creating a subset of wcvp where we only have the accepted species in the order
cat("Creating a subset of wcvp where we only have the accepted species in the order \n")
wcvp_accepted_species_orders_subset <- wcvp_accepted_species_orders[which(wcvp_accepted_species_orders$order %in% order_in_question),]

# Now that i have a wcvp with added orders, I should be able to use this to slice the TACTed tree into each order
tree_order <- keep.tip(tree, tree$tip.label[which(tree$tip.label %in% wcvp_accepted_species_orders_subset$plant_name_id)])
length(tree_order$tip.label)

tree$tip.label[which(tree$tip.label %in% wcvp_accepted_species_orders_subset$plant_name_id)]

head(wcvp_accepted_species_orders_subset)

all(tree$tip.label %in% wcvp_accepted_species$plant_name_id)
length(which(tree$tip.label %in% wcvp_accepted_species$plant_name_id))
length(tree$tip.label)
